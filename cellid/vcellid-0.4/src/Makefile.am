## Process this file with automake to produce Makefile.in

INCLUDES = \
	-DPACKAGE_DATA_DIR=\""$(datadir)"\" \
	-DPACKAGE_LOCALE_DIR=\""$(prefix)/$(DATADIRNAME)/locale"\" \
	@PACKAGE_CFLAGS@

bin_PROGRAMS = vcellid

vcellid_SOURCES = \
	main.c main.h\
	support.c support.h \
	interface.c interface.h \
	callbacks.c callbacks.h

vcellid_LDADD = @PACKAGE_LIBS@ $(INTLLIBS) -lgtkimageview

