AUTOMAKE_OPTIONS=foreign
EXTRA_DIST = $(top_srcdir)/deprecated/pymoab/* $(srcdir)/* $(top_srcdir)/test/pymoab/*
EXTRA_DIST += $(srcdir)/interface_notes.txt
test_files = test_core.py test_hom.py test_scd.py test_range.py test_tag.py test_topo_util.py test_skinner.py
local_install_path = $(subst $(prefix)/,,$(pythondir))
all-local: local-install
local-install:
	$(MKDIR_P) $(top_builddir)/$(local_install_path)
	(CFLAGS='-w -fPIC' PYTHONPATH=$(top_builddir)/$(local_install_path)/:$(PYTHONPATH) \
		$(PYTHON) -m pip install -e $(top_builddir) --prefix $(top_builddir)/$(local_install_path))
install-exec-local: all-local
	$(MKDIR_P) $(DESTDIR)$(pythondir)
	(PYTHONPATH=$(DESTDIR)$(pythondir):$(PYTHONPATH) $(PYTHON) -m pip install $(top_builddir) \
        --prefix $(DESTDIR)$(prefix) )
check: all-local
	cd $(top_srcdir)/test/pymoab; \
	for i in $(test_files); \
	do \
		(@SLIBENVPATH@=$(top_builddir)/src/.libs/:$(@SLIBENVPATH@) PYTHONPATH=$(top_builddir):$(PYTHONPATH) $(PYTHON) $$i) || exit 1; \
	done
clean-local:
	-rm -rf $(top_srcdir)/build $(top_srcdir)/lib *.cpp __pycache__ *.so *.pyc $(top_srcdir)/test/pymoab/*.pyc $(top_srcdir)/pymoab.egg-info
uninstall:
	rm -rf $(DESTDIR)$(pythondir)/*
CLEANFILES = $(top_srcdir)/test/pymoab/quads_test.vtk $(top_srcdir)/test/pymoab/outfile.h5m $(top_srcdir)/test/pymoab/write_tag_test.h5m $(local_install_path)/*