
SUBDIRS =

noinst_PROGRAMS =  stationary 

stationary_SOURCES = stationary.cc

stationary_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(UG_CPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	$(SUPERLU_CPPFLAGS)
# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
stationary_LDADD = \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(SUPERLU_LIBS) \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	$(UG_LDFLAGS) $(UG_LIBS) \
	$(DUNEMPILIBS)	\
	$(LDADD) \
	-lsuperlu
stationary_LDFLAGS = $(AM_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(UG_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	$(SUPERLU_LDFLAGS) \
	$(DUNE_LDFLAGS)



# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules

