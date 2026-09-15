# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include makefile.config

default: build

# advanced/smoothing is listed separately rather than walked from
# advanced/makefile, which has no SUBDIRS of its own; until now nothing
# descended into it at all.
SUBDIRS = beginner basic intermediate advanced advanced/smoothing fortran earthsystem# cmake

build: ${MOAB_LIBFILE}
	@for exampledir in ${SUBDIRS} ; do \
		${MAKE} -C $$exampledir all ; \
	done

run:
	@for exampledir in ${SUBDIRS} ; do \
		${MAKE} -C $$exampledir run ; \
	done

clean: clobber
	@for exampledir in ${SUBDIRS} ; do \
		${MAKE} -C $$exampledir clean ; \
	done

