# Automake template used in the autoconf/automake build system
# See man/info pages for autoconf and automake
# Also, see man/info pages for pkg-config

# Add Fortran .mod files to the list of things to be cleaned

CLEANFILES = *.mod

# Put the cosmological survey model code into library

lib_LIBRARIES = libcosmo.a

libcosmo_a_SOURCES = \
	quadpack.f90 \
	quadrature.f90 \
	utils.f90 \
	constants.f90 \
	cosmology.f90 \
	structure.f90

# Build some executables that use the library

bin_PROGRAMS = cosmosurvey cosmoplot cosmotest

# cosmosurvey
cosmosurvey_SOURCES = \
	kraken.f90  \
	survey2.f90 \
	survey3.f90 \
	cosmosurvey.f90
cosmosurvey_LDADD = libcosmo.a

# cosmoplot
cosmoplot_FCFLAGS = $(PLPLOT_CFLAGS)
cosmoplot_SOURCES = \
	kraken.f90  \
	survey2.f90 \
	survey3.f90 \
	cosmoplot.f90
cosmoplot_LDADD = libcosmo.a $(PLPLOT_LIBS)

# cosmotest
cosmotest_FCFLAGS = $(PLPLOT_CFLAGS)
cosmotest_SOURCES = \
	kraken.f90  \
	survey2.f90 \
	survey3.f90 \
	cosmotest.f90
cosmotest_LDADD = libcosmo.a $(PLPLOT_LIBS)

