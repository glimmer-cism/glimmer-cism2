SUBDIRS = m4macros doc scripts libglimmer libglimmer-solve libglide libglint \
          utils tests example-drivers extensions
if BUILD_PYTHON
  SUBDIRS += PyGC
endif
EXTRA_DIST = make.rules glimmer-cism-config.in
bin_SCRIPTS=glimmer-cism-config

datadir = tests example-drivers

data::
	list='$(datadir)'; for dir in $$list; do \
          (cd $$dir && $(MAKE) $@); \
        done
