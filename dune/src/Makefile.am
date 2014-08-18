
SUBDIRS = data Bar1AlphaStationary Bar1AlphaInstationary Elston-Model 2d3d-comparison

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-fem=$(DUNE_FEM_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules
