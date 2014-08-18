# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src])
AM_CONFIG_HEADER([config.h])


# we need no more than the standard DE-stuff
# this module depends on dune-common dune-grid dune-istl dune-localfunctions dune-pdelab
# this implies checking for [dune-common], [dune-grid], [dune-istl], [dune-localfunctions], [dune-pdelab]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  src/data/Makefile
  src/Elston-Model/Makefile
  src/Elston-Model/bc/Makefile	
  src/Bar1AlphaStationary/Makefile
  src/Bar1AlphaStationary/bc/Makefile	
  src/Bar1AlphaInstationary/Makefile
  src/Bar1AlphaInstationary/bc/Makefile
  src/2d3d-comparison/Makefile
  src/2d3d-comparison/bc/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/ExtracellularCommunication/Makefile
  m4/Makefile
  ExtracellularCommunication.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
