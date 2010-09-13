bin_PROGRAMS  =  simple_glide
if BUILD_EXT_EROSION
 bin_PROGRAMS += simple_erosion
endif

lib_LTLIBRARIES = libglimmer-simple.la

libglimmer_simple_la_SOURCES = simple_forcing.F90
libglimmer_simple_la_LDFLAGS=

simple_glide_SOURCES = simple_glide.F90
if USE_GLIDE_IF
simple_glide_LDADD = libglimmer-simple.la $(top_builddir)/libglide_if/libglide_if.la $(top_builddir)/libglide/libglide.la $(top_builddir)/libglimmer-solve/libglimmer-solve.la \
                    $(top_builddir)/libglimmer/libglimmer-IO.la $(top_builddir)/libglimmer/libglimmer.la \
                    $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS) $(BLAS_LIBS)
else
simple_glide_LDADD = libglimmer-simple.la $(top_builddir)/libglide/libglide.la $(top_builddir)/libglimmer-solve/libglimmer-solve.la \
                    $(top_builddir)/libglimmer/libglimmer-IO.la $(top_builddir)/libglimmer/libglimmer.la \
                    $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS) $(BLAS_LIBS)
endif
simple_erosion_SOURCES = simple_erosion.F90
if USE_GLIDE_IF
simple_erosion_LDADD = libglimmer-simple.la  $(top_builddir)/extensions/erosion/src/liberosion.la $(top_builddir)/libglide_if/libglide_if.la $(top_builddir)/libglide/libglide.la \
                    $(top_builddir)/libglimmer-solve/libglimmer-solve.la $(top_builddir)/libglimmer/libglimmer-IO.la $(top_builddir)/libglimmer/libglimmer.la \
                    $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS) $(BLAS_LIBS)
else
simple_erosion_LDADD = libglimmer-simple.la  $(top_builddir)/extensions/erosion/src/liberosion.la $(top_builddir)/libglide/libglide.la \
                    $(top_builddir)/libglimmer-solve/libglimmer-solve.la $(top_builddir)/libglimmer/libglimmer-IO.la $(top_builddir)/libglimmer/libglimmer.la \
                    $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS) $(BLAS_LIBS)
endif

AM_FCFLAGS = $(NETCDF_FCFLAGS) $(F90_MODFLAG)$(top_builddir)/libglimmer-solve $(F90_MODFLAG)$(top_builddir)/libglide $(F90_MODFLAG)$(top_builddir)/libglimmer \
             $(F90_MODFLAG)$(top_builddir)/extensions/erosion/src/
if USE_GLIDE_IF
AM_FCFLAGS += $(F90_MODFLAG)$(top_builddir)/libglide_if
endif


CLEANFILES = f90_dep.mak *.mod
-include f90_dep.mak
include $(top_srcdir)/extra_rules.am
