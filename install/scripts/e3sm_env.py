#!/usr/bin/env python3
"""
Emit a sourceable bash snippet from E3SM's config_machines.xml for a given
(machine, compiler) pair.

Used by install-moab.sh under --profile=e3sm to load the same modules and
environment variables that an E3SM build on the target machine would use,
without duplicating that knowledge inside install-moab.sh.

Output goes to stdout. Diagnostics go to stderr. Exit codes:
    0  success
    2  argument / lookup failure (machine not in XML, --e3sm-root missing, ...)
    3  XML parse / unexpected error

Two-tier strategy:
    Tier 1: prefer CIME's Machines class for XML traversal (handles schema
            evolution and exposes machine_node / get_children() helpers).
    Tier 2: fall back to xml.etree.ElementTree if CIME isn't importable
            (refactor in a future E3SM, or no cime/ dir).

Variable substitution: we deliberately do NOT use CIME's get_resolved_value()
because it eagerly resolves $ENV{}/$SHELL{} against the LOCAL env (the host
generating the script), which is wrong when generating on a login node /
laptop and consuming on a different machine. Instead we translate:
    $ENV{VAR}     -> ${VAR:-}      (safe under set -u; resolves at source time)
    $SHELL{cmd}   -> $(cmd)        (resolves at source time on the target)
This keeps the emitted snippet portable and idempotent.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable, Optional


# ---------- variable substitution -------------------------------------------------

_RE_SHELL = re.compile(r"\$SHELL\{([^}]*)\}")     # nested {} not used in any current entry
_RE_ENV = re.compile(r"\$ENV\{([A-Za-z_][A-Za-z0-9_]*)\}")


def translate_value(raw: str) -> str:
    """CIME-syntax -> bash-syntax. Leaves unrelated $VAR and `$()` alone."""
    if raw is None:
        return ""
    out = _RE_SHELL.sub(lambda m: f"$({m.group(1)})", raw)
    out = _RE_ENV.sub(lambda m: "${" + m.group(1) + ":-}", out)
    # Escape backslashes and double quotes for safe inclusion inside `export VAR="..."`
    out = out.replace("\\", "\\\\").replace('"', r"\"")
    return out


def block_matches(attrs: dict, compiler: str, mpilib: str) -> bool:
    """A compiler/mpilib-attributed block applies if its attrs match the request.
    Blocks with BUILD_THREADED are skipped (we don't build MOAB threaded here)."""
    if "BUILD_THREADED" in attrs:
        return False
    if "compiler" in attrs and attrs["compiler"] != compiler:
        return False
    if "mpilib" in attrs and attrs["mpilib"] != mpilib:
        return False
    # debug_compiler, build_threaded etc. handled implicitly by the attr filter
    return True


# ---------- Tier 1: CIME API ------------------------------------------------------


def emit_via_cime(e3sm_root: Path, machine: str, compiler: str, mpilib: str,
                  shell_lang: str = "sh") -> None:
    """Walk the XML via CIME's Machines class. Raises on lookup failure."""
    sys.path.insert(0, str(e3sm_root / "cime"))
    from CIME.XML.machines import Machines  # type: ignore

    m = Machines(machine=machine)
    mach_node = m.machine_node

    # Resolve mpilib if not specified: pick first listed in <MPILIBS>
    if not mpilib:
        mpilibs = m.get_value("MPILIBS") or ""
        mpilib = mpilibs.split(",")[0].strip() if mpilibs else ""

    # Validate compiler is supported on this machine (warning only -- user may know better)
    compilers = (m.get_value("COMPILERS") or "").split(",")
    if compiler not in [c.strip() for c in compilers if c.strip()]:
        sys.stderr.write(
            f"# WARN: compiler '{compiler}' not in {machine}'s COMPILERS={compilers}; proceeding\n"
        )

    print(f"# === e3sm_env.py output: machine={machine} compiler={compiler} mpilib={mpilib} ===")
    print(f"# (via CIME API at {e3sm_root}/cime)")
    print()

    # Module system: source the init script, then issue commands
    ms_node = m.get_child("module_system")
    if ms_node is not None:
        init_node = m.get_optional_child("init_path", attributes={"lang": shell_lang}, root=ms_node)
        cmd_node = m.get_optional_child("cmd_path", attributes={"lang": shell_lang}, root=ms_node)
        init_path = m.text(init_node) if init_node is not None else ""
        cmd_path = m.text(cmd_node) if cmd_node is not None else "module"
        if init_path:
            print(f"# Lmod / module-system init")
            print(f"if [ -f '{init_path}' ]; then")
            print(f"    . '{init_path}'")
            print(f"else")
            print(f"    echo 'WARN [e3sm_env]: module init script not found at {init_path}' >&2")
            print(f"fi")
            print()

        # Modules
        emitted_any = False
        for blk in m.get_children("modules", root=ms_node):
            attrs = m.attrib(blk) or {}
            if not block_matches(attrs, compiler, mpilib):
                continue
            if not emitted_any:
                print("# Module loads")
                emitted_any = True
            for cmd in m.get_children(root=blk):
                action = (m.attrib(cmd) or {}).get("name", "load")
                args_str = (m.text(cmd) or "").strip()
                if not args_str:
                    continue
                # `module purge` etc. take no args
                if action in ("purge",):
                    print(f"{cmd_path} {action}")
                elif action == "swap":
                    # e.g. <command name="swap">PrgEnv-pgi PrgEnv-gnu</command>
                    print(f"{cmd_path} swap {args_str}")
                else:
                    print(f"{cmd_path} {action} {args_str}")
        if emitted_any:
            print()

    # Environment variables
    emitted_env_header = False
    saw_hdf5_root = False
    for blk in m.get_children("environment_variables", root=mach_node):
        attrs = m.attrib(blk) or {}
        if not block_matches(attrs, compiler, mpilib):
            continue
        if not emitted_env_header:
            print("# Environment variables")
            emitted_env_header = True
        for env in m.get_children("env", root=blk):
            name = (m.attrib(env) or {}).get("name")
            raw = m.text(env)
            if not name:
                continue
            if name == "HDF5_ROOT":
                saw_hdf5_root = True
            translated = translate_value(raw)
            print(f'export {name}="{translated}"')
    if emitted_env_header:
        print()

    # HDF5_ROOT gap diagnostic. Bebop / Improv / Cray machines deliberately
    # don't export HDF5_ROOT in config_machines.xml, but install-moab.sh
    # requires it. Surface this clearly so the user knows to pass --hdf5-root=.
    if not saw_hdf5_root:
        print(f"# NOTE [e3sm_env]: machine '{machine}' does not export HDF5_ROOT in config_machines.xml.")
        print(f"# If install-moab.sh fails with 'HDF5_ROOT is not set', either:")
        print(f"#   * pass --hdf5-root=PATH explicitly, or")
        print(f"#   * load a separate hdf5 module that exports it, or")
        print(f"#   * for Cray PrgEnv: HDF5 lives at \\$CRAY_HDF5_PARALLEL_PREFIX after `module load cray-hdf5-parallel`")
        print()
    print(f"# === end e3sm_env.py output ===")


# ---------- Tier 2: direct XML fallback -------------------------------------------


def emit_via_xml(e3sm_root: Path, machine: str, compiler: str, mpilib: str,
                 shell_lang: str = "sh") -> None:
    """Pure-stdlib fallback: walk config_machines.xml with ElementTree."""
    import xml.etree.ElementTree as ET

    xml_path = e3sm_root / "cime_config" / "machines" / "config_machines.xml"
    tree = ET.parse(xml_path)
    root = tree.getroot()
    mach_node = root.find(f"./machine[@MACH='{machine}']")
    if mach_node is None:
        sys.stderr.write(f"ERROR: machine '{machine}' not found in {xml_path}\n")
        sys.exit(2)

    if not mpilib:
        mpilibs = (mach_node.findtext("MPILIBS") or "").strip()
        mpilib = mpilibs.split(",")[0].strip() if mpilibs else ""

    compilers = [c.strip() for c in (mach_node.findtext("COMPILERS") or "").split(",") if c.strip()]
    if compiler not in compilers:
        sys.stderr.write(
            f"# WARN: compiler '{compiler}' not in {machine}'s COMPILERS={compilers}; proceeding\n"
        )

    print(f"# === e3sm_env.py output: machine={machine} compiler={compiler} mpilib={mpilib} ===")
    print(f"# (via direct XML fallback; CIME API was not available)")
    print()

    ms_node = mach_node.find("module_system")
    if ms_node is not None:
        init_node = ms_node.find(f"./init_path[@lang='{shell_lang}']")
        cmd_node = ms_node.find(f"./cmd_path[@lang='{shell_lang}']")
        init_path = (init_node.text if init_node is not None else "") or ""
        cmd_path = (cmd_node.text if cmd_node is not None else "module") or "module"
        if init_path:
            print(f"# Lmod / module-system init")
            print(f"if [ -f '{init_path}' ]; then")
            print(f"    . '{init_path}'")
            print(f"else")
            print(f"    echo 'WARN [e3sm_env]: module init script not found at {init_path}' >&2")
            print(f"fi")
            print()

        emitted_any = False
        for blk in ms_node.findall("modules"):
            if not block_matches(blk.attrib or {}, compiler, mpilib):
                continue
            if not emitted_any:
                print("# Module loads")
                emitted_any = True
            for cmd in list(blk):
                action = cmd.attrib.get("name", "load")
                args_str = (cmd.text or "").strip()
                if not args_str and action != "purge":
                    continue
                if action == "purge":
                    print(f"{cmd_path} {action}")
                elif action == "swap":
                    print(f"{cmd_path} swap {args_str}")
                else:
                    print(f"{cmd_path} {action} {args_str}")
        if emitted_any:
            print()

    emitted_env_header = False
    saw_hdf5_root = False
    for blk in mach_node.findall("environment_variables"):
        if not block_matches(blk.attrib or {}, compiler, mpilib):
            continue
        if not emitted_env_header:
            print("# Environment variables")
            emitted_env_header = True
        for env in blk.findall("env"):
            name = env.attrib.get("name")
            raw = env.text
            if not name:
                continue
            if name == "HDF5_ROOT":
                saw_hdf5_root = True
            translated = translate_value(raw)
            print(f'export {name}="{translated}"')
    if emitted_env_header:
        print()

    if not saw_hdf5_root:
        print(f"# NOTE [e3sm_env]: machine '{machine}' does not export HDF5_ROOT in config_machines.xml.")
        print(f"# If install-moab.sh fails with 'HDF5_ROOT is not set', either:")
        print(f"#   * pass --hdf5-root=PATH explicitly, or")
        print(f"#   * load a separate hdf5 module that exports it, or")
        print(f"#   * for Cray PrgEnv: HDF5 lives at \\$CRAY_HDF5_PARALLEL_PREFIX after `module load cray-hdf5-parallel`")
        print()
    print(f"# === end e3sm_env.py output ===")


# ---------- main ------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--e3sm-root", required=True, type=Path,
                    help="Path to E3SM checkout (must contain cime_config/machines/config_machines.xml)")
    ap.add_argument("--machine", required=True,
                    help="Machine name as it appears in config_machines.xml (e.g. bebop, pm-cpu)")
    ap.add_argument("--compiler", required=True,
                    help="Compiler family (gnu, intel, cray, nvidia, aocc, ...)")
    ap.add_argument("--mpilib", default="",
                    help="MPI library (auto-detected from MPILIBS if omitted)")
    ap.add_argument("--shell-lang", default="sh", choices=["sh", "csh", "python"],
                    help="Shell flavour for module init script (default: sh)")
    ap.add_argument("--no-cime", action="store_true",
                    help="Skip the CIME API tier; use direct XML parse only (for testing the fallback)")
    args = ap.parse_args()

    e3sm_root = args.e3sm_root.resolve()
    cm_xml = e3sm_root / "cime_config" / "machines" / "config_machines.xml"
    if not cm_xml.is_file():
        sys.stderr.write(f"ERROR: not an E3SM checkout: {cm_xml} does not exist\n")
        return 2

    if not args.no_cime:
        try:
            emit_via_cime(e3sm_root, args.machine, args.compiler, args.mpilib,
                          shell_lang=args.shell_lang)
            return 0
        except SystemExit:
            raise
        except ImportError as e:
            sys.stderr.write(f"# Tier 1 (CIME API) unavailable ({e}); falling back to direct XML parse\n")
        except Exception as e:
            sys.stderr.write(
                f"# Tier 1 (CIME API) failed ({type(e).__name__}: {e}); falling back to direct XML parse\n"
            )

    try:
        emit_via_xml(e3sm_root, args.machine, args.compiler, args.mpilib,
                     shell_lang=args.shell_lang)
        return 0
    except SystemExit:
        raise
    except Exception as e:
        sys.stderr.write(f"ERROR [e3sm_env, XML fallback]: {type(e).__name__}: {e}\n")
        return 3


if __name__ == "__main__":
    sys.exit(main())
