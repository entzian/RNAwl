
bin_PROGRAMS = RNAwl

SUFFIXES = _cmdl.c _cmdl.h .ggo
GENGETOPT_CMDL =  wanglandau_cmdl.c wanglandau_cmdl.h
GENGETOPT_FILES =  wanglandau.ggo
EXTRA_DIST =  ${GENGETOPT_FILES} ${GENGETOPT_CMDL} 
# remove all generated content
maintainer-clean-local:
	rm -rf $(GENGETOPT_CMDL)

# how we make the commandline interface parser
.ggo_cmdl.c:
	$(AM_V_GEN)gengetopt -i $<

_cmdl.c_cmdl.h :

RNAwl_SOURCES =	${GENGETOPT_CMDL} \
            main.c\
			moves.c\
			wl_options.c\
			wanglandau.c\
			wl_rna.c
			

AM_CFLAGS = ${GSL_CFLAGS} ${ViennaRNA_CFLAGS} -g3 -O0
AM_CPPFLAGS = -I${includedir} -I.

LDADD = ${GSL_LIBS} ${ViennaRNA_LIBS}

