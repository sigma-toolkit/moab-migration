# mbtagprop {#tools_mbtagprop}

Propagate tags from entity sets to the mesh entities those sets contain. Common use case: lift a material/boundary identifier from a `MATERIAL_SET` set onto every element it owns.

## Synopsis

```
mbtagprop <options> <input_file> <output_file>
```

## Options

| Flag | Description |
|---|---|
| `-t <ident_tag>[=<value>]` | Identifier tag used to select entity sets. **At least one is required.** May be repeated; logical OR over the matches. |
| `-d <data_tag>[=<default>]` | Existing tag whose value is propagated onto contained entities. If omitted, the value of `ident_tag` is propagated. The optional `<default>` is applied to entities whose owning set has no `data_tag`. |
| `-c <data_tag=type:size>[=<default>]` | Same as `-d`, but **creates** the `data_tag` rather than requiring it to exist. |
| `-w <write_tag>` | Rename the propagated tag in the output. Defaults to `data_tag`. |
| `-n` | Propagate to **nodes** in the identified sets only. |
| `-e` | Propagate to **elements** in the identified sets only. |
| *(neither `-n` nor `-e`)* | Propagate to both nodes and elements *(default)*. |

If multiple `-t` flags are given, `-d` (or `-c`) **must** be specified explicitly — there's no implicit "use the ident tag" fallback when the ident tag is ambiguous.

## Tag value syntax

Tags are specified as `<name>=[value]`; the value is optional.

### Value formats

- **Opaque tags**: a numeric value may be specified in hexadecimal with a `0x` prefix. Without `0x` the value is treated as a string and NULL-padded to the tag size.
- **Integral tags** (`INTEGER`, `BIT`, `HANDLE`): parsed as integers; standard `0x` (hex) and `0` (octal) prefixes are accepted.
- **`DOUBLE` tags**: base-10 only; exponential notation (`1.5e-3`) is accepted.
- **Array-valued tags** (non-opaque): comma-separated list, no spaces. E.g. `tag=1,2,3`.

### Creating a new tag

`-c` requires a type and size as part of the spec:

```
-c <name>=<type>:<size>[=<default_value>]
```

- `<type>` ∈ `int | double | opaque | handle | bit`
- `<size>` = number of values of `<type>` (or number of bytes for `opaque`)
- `<default_value>` optional, follows the same value-syntax rules above

## Examples

```bash
# Lift MATERIAL_SET tag values onto all owned entities
mbtagprop -t MATERIAL_SET  in.h5m  out.h5m

# Lift only onto the elements (not nodes) of NEUMANN_SET == 7
mbtagprop -t NEUMANN_SET=7 -e  in.h5m  out.h5m

# Create a new integer tag "block_id" and propagate from MATERIAL_SET values
mbtagprop -t MATERIAL_SET -c block_id=int:1  in.h5m  out.h5m

# Multiple ident tags + explicit data tag
mbtagprop -t TAG1 -t TAG2 -d shared_tag  in.h5m  out.h5m
```
