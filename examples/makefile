# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include makefile.config

default: build

SUBDIRS = beginner basic intermediate advanced fortran earthsystem# cmake

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

