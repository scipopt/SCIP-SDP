#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*                                                                           */
#/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
#/* semidefinite programms based on SCIP.                                     */
#/*                                                                           */
#/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
#/*                         EDOM, FAU Erlangen-Nürnberg                       */
#/*               2014-2016 Discrete Optimization, TU Darmstadt               */
#/*                                                                           */
#/*                                                                           */
#/* This program is free software; you can redistribute it and/or             */
#/* modify it under the terms of the GNU Lesser General Public License        */
#/* as published by the Free Software Foundation; either version 3            */
#/* of the License, or (at your option) any later version.                    */
#/*                                                                           */
#/* This program is distributed in the hope that it will be useful,           */
#/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
#/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
#/* GNU Lesser General Public License for more details.                       */
#/*                                                                           */
#/* You should have received a copy of the GNU Lesser General Public License  */
#/* along with this program; if not, write to the Free Software               */
#/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
#/*                                                                           */
#/*                                                                           */
#/* Based on SCIP - Solving Constraint Integer Programs                       */
#/* Copyright (C) 2002-2016 Zuse Institute Berlin                             */
#/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
#/* see file COPYING in the SCIP distribution.                                */
#/*                                                                           */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#@file    Makefile
#@brief   Makefile for C++ SDP-Interface for SCIP
#@author  Sonja Mars, Lars Schewe, Marc Pfetsch, Tristan Gally, Ambros Gleixner

#-----------------------------------------------------------------------------
# own variables
#-----------------------------------------------------------------------------

SCIPSDPVERSION	=	2.1.0
SDPS		=	none

GCCWARN		+= 	-Wextra

# parameters
TIME     	=  	3600
NODES           =       2100000000
MEM		=	6144
THREADS         =       1
PERMUTE         =       0
DISPFREQ	=	10000
FEASTOL		=	default
TEST		=	short
SETTINGS        =       default
CONTINUE	=	false
LOCK		=	false
VALGRIND	=	false
CLIENTTMPDIR    =       /tmp
OPTCOMMAND	=	optimize
SOFTLINKS	=
MAKESOFTLINKS	=	true
OPENBLAS	=	true


#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------

# save directory to be able to locate library files
ifeq ($(OSTYPE),mingw)
SCIPSDPDIR	=	./
else
SCIPSDPDIR	=	$(realpath .)
endif
SCIPSDPLIBDIR	=	lib

SCIPDIR		= 	$(SCIPSDPDIR)/lib/scip

# check whether SCIPDIR exists
ifeq ("$(wildcard $(SCIPDIR))","")
$(error Please add a soft-link to the SCIP directory as $(SCIPSDPDIR)/lib/scip)
endif

# load SCIP project makefile
include $(SCIPDIR)/make/make.project


#-----------------------------------------------------------------------------
# settings for SDP solver
#-----------------------------------------------------------------------------

LDFLAGS 	+= 	-lobjscip

SDPIOPTIONS	=
SDPIINC		=
SDPILIB		=
SDPICSRC	=
SDPICCSRC	=
SDPICCOBJ	=


#-----------------------------------------------------------------------------
# DSDP solver
SDPIOPTIONS	+=	dsdp
ifeq ($(SDPS),dsdp)
SDPILIB		= 	-L$(SCIPSDPLIBDIR) -ldsdp -llapack -lblas
SDPIINC		= 	-I$(SCIPSDPLIBDIR)/dsdpinc
SDPICSRC 	= 	src/sdpi/sdpisolver_dsdp.c \
			src/sdpi/lapack_dsdp.c
SDPIOBJ 	= 	$(OBJDIR)/sdpi/sdpisolver_dsdp.o \
			$(OBJDIR)/sdpi/lapack_dsdp.o
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/dsdpinc
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/libdsdp.$(STATICLIBEXT)
SDPIINSTMSG	=	"  -> \"dsdpinc\" is the path to the DSDP \"include\" directory, e.g., \"<DSDP-path>/include\".\n"
SDPIINSTMSG	+=	" -> \"libdsdp.*\" is the path to the DSDP library, e.g., \"<DSDP-path>/lib/libdsdp.$(STATICLIBEXT)\""
endif

#-----------------------------------------------------------------------------
# SDPA solver
SDPIOPTIONS	+=	sdpa
ifeq ($(SDPS),sdpa)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/sdpainc
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/libsdpa.$(STATICLIBEXT)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/libsdpa.$(SHAREDLIBEXT)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/mumpsinc
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/mumpslib
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/mumpslibseq
ifeq ($(OPENBLAS),true)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/libopenblas.$(SHAREDLIBEXT).0
endif
SDPIINSTMSG	=	"  -> \"sdpainc\" is the path to the SDPA \"include\" directory, e.g., \"<SDPA-path>/include\".\n"
SDPIINSTMSG	+=	" -> \"libsdpa.*\" is the path to the SDPA library, e.g., \"<SDPA-path>/lib/libsdpa.a\".\n"
SDPIINSTMSG	+=	" -> \"mumpsinc\" is the path to the Mumps \"include\" directory, e.g., \"<SDPA-path>/mumps/build/include\".\n"
SDPIINSTMSG	+=	" -> \"mumpslib\" is the path to the Mumps \"lib\" directory, e.g., \"<SDPA-path>/mumps/build/lib\".\n"
SDPIINSTMSG	+=	" -> \"mumpslibseq\" is the path to the Mumps \"libseq\" directory, e.g., \"<SDPA-path>/mumps/build/libseq\".\n"
ifeq ($(OPENBLAS),true)
SDPIINSTMSG	+=	" -> \"libopenblas.$(SHAREDLIBEXT).0\" is the openblas library.\n"
SDPILIB		=      	-L$(SCIPSDPLIBDIR) -lsdpa -L$(SCIPSDPLIBDIR)/mumpslib -ldmumps -lmumps_common -lpord -L$(SCIPSDPLIBDIR)/mumpslibseq -lmpiseq \
			$(SCIPSDPLIBDIR)/libopenblas.$(SHAREDLIBEXT).0 -Wl,-rpath,$(SCIPSDPDIR)/$(SCIPSDPLIBDIR) \
			-lgfortran -L/lib/x86_64-linux-gnu -lpthread -lgomp
else
SDPILIB		=      	-L$(SCIPSDPLIBDIR) -lsdpa -L$(SCIPSDPLIBDIR)/mumpslip -ldmumps -lmumps_common -lpord -L$(SCIPSDPLIBDIR)/mumpslibseq -lmpiseq \
			-lgfortran -llapack -lblas
endif
SDPIINC		=      	-I$(SCIPSDPLIBDIR)/sdpainc -I$(SCIPSDPLIBDIR)/mumpsinc
SDPICCSRC 	= 	src/sdpi/sdpisolver_sdpa.cpp
SDPICSRC	=	src/sdpi/lapack_sdpa.c
SDPIOBJ 	= 	$(OBJDIR)/sdpi/sdpisolver_sdpa.o $(OBJDIR)/sdpi/lapack_sdpa.o
endif

#-----------------------------------------------------------------------------
# no solver
SDPIOPTIONS	+=	none
ifeq ($(SDPS),none)
SDPILIB		= 	-L$(SCIPSDPLIBDIR) -ldsdp -llapack -lblas
SDPICSRC 	= 	src/sdpi/sdpisolver_none.c src/sdpi/lapack_dsdp.c
SDPIOBJ 	= 	$(OBJDIR)/sdpi/sdpisolver_none.o $(OBJDIR)/sdpi/lapack_dsdp.o
SETTINGS	= 	lp_approx
endif

# include install/uninstall targets
-include make/make.install

LINKSMARKERFILE	=	$(LIBDIR)/linkscreated.$(LPS)-$(LPSOPT).$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX).$(ZIMPL)-$(ZIMPLOPT).$(IPOPT)-$(IPOPTOPT).$(GAMS)

#-----------------------------------------------------------------------------

SDPOBJSUBDIRS	=	$(OBJDIR)/scipsdp \
			$(OBJDIR)/sdpi


#-----------------------------------------------------------------------------
# main program
#-----------------------------------------------------------------------------

MAINNAME	=	scipsdp
MAINCOBJ	=	scipsdp/SdpVarmapper.o \
			scipsdp/SdpVarfixer.o \
			scipsdp/cons_sdp.o \
			scipsdp/cons_savedsdpsettings.o \
			scipsdp/relax_sdp.o \
			scipsdp/disp_sdpiterations.o \
			scipsdp/disp_sdpavgiterations.o \
			scipsdp/disp_sdpfastsettings.o \
			scipsdp/disp_sdppenalty.o \
			scipsdp/disp_sdpunsolved.o \
			scipsdp/prop_sdpredcost.o \
			scipsdp/branch_sdpmostfrac.o \
			scipsdp/branch_sdpmostinf.o \
			scipsdp/branch_sdpobjective.o \
			scipsdp/branch_sdpinfobjective.o \
			scipsdp/heur_sdpfracdiving.o \
			scipsdp/heur_sdprand.o \
			scipsdp/prop_sdpobbt.o \
			sdpi/sdpi.o

MAINCCOBJ 	=	scipsdp/main.o \
			scipsdp/objreader_sdpa.o \
			scipsdp/ScipStreamBuffer.o


MAINCSRC	=	$(addprefix $(SRCDIR)/,$(MAINCOBJ:.o=.c))
MAINCCSRC 	=	$(addprefix $(SRCDIR)/,$(MAINCCOBJ:.o=.cpp))
MAINDEP 	=	$(SRCDIR)/depend.cppmain.$(OPT)

# @todo possibly add LPS
MAINFILE	=	$(BINDIR)/$(MAINNAME).$(BASE).$(SDPS)$(EXEEXTENSION)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINCOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINCOBJ))
MAINCCOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINCCOBJ))

ALLSRC		=	$(MAINCSRC) $(MAINCCSRC) $(SDPICSRC) $(SDPICCSRC)
LINKSMARKERFILE =	$(SCIPSDPLIBDIR)/linkscreated.$(SDPS).$(LPS)-$(LPSOPT).$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
LASTSETTINGS 	=	$(OBJDIR)/make.lastsettings


#-----------------------------------------------------------------------------
# rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINCOBJFILES) $(MAINCCOBJFILES) $(SDPIOBJ) $(MAINSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: checkdefines
checkdefines:
ifeq ($(SDPIOBJ),)
		$(error invalid SDP solver selected: SDPS=$(SDPIS). Possible options are: $(SDPIOPTIONS))
endif

.PHONY: preprocess
preprocess:     checkdefines
		@$(SHELL) -ec 'if test ! -e $(LINKSMARKERFILE) ; \
			then \
				echo "-> generating necessary links" ; \
				$(MAKE) -j1 $(LINKSMARKERFILE) ; \
			fi'
		@$(MAKE) touchexternal

.PHONY: tags
tags:
		rm -f TAGS; ctags -e src/*/*.c src/*/*.cpp src/*/*.h $(SCIPDIR)/src/*/*.c $(SCIPDIR)/src/*/*.h;

.PHONY: lint
lint:		$(MAINCSRC) $(MAINCCSRC) $(SDPICSRC) $(SDPICCSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) -i lint co-gcc.lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'

.PHONY: doc
doc:
		cd doc; $(DOXY) $(MAINNAME).dxy

$(MAINSHORTLINK): $(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@mkdir -p $(OBJDIR);

$(SDPOBJSUBDIRS):	| $(OBJDIR)
		@-mkdir -p $(SDPOBJSUBDIRS);

$(SCIPSDPLIBDIR):
		@-mkdir -p $(SCIPSDPLIBDIR)

$(BINDIR):
		-@test -d $(BINDIR) || { \
		echo "-> Creating $(BINDIR) directory"; \
		mkdir -p $(BINDIR); }

.PHONY: clean
clean:
ifneq ($(OBJDIR),)
		@-rm -f $(LASTSETTINGS)
		@-rm -f $(OBJDIR)/scipsdp/*.o
		@-rm -f $(OBJDIR)/sdpi/*.o
		@-rmdir $(OBJDIR)/scipsdp
	 	@-rmdir $(OBJDIR)/sdpi
		@-rmdir $(OBJDIR)
endif
		-rm -f $(MAINFILE)

#-----------------------------------------------------------------------------
-include $(LASTSETTINGS)

.PHONY: touchexternal
touchexternal:	$(SCIPSDPLIBDIR) $(OBJDIR)
ifneq ($(SHARED),$(LAST_SHARED))
		@-touch $(ALLSRC)
endif
ifneq ($(USRFLAGS),$(LAST_USRFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USROFLAGS),$(LAST_USROFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USRCFLAGS),$(LAST_USRCFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USRCXXFLAGS),$(LAST_USRCXXFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USRLDFLAGS),$(LAST_USRLDFLAGS))
		@-touch -c $(ALLSRC)
endif
ifneq ($(USRARFLAGS),$(LAST_USRARFLAGS))
		@-touch -c $(ALLSRC)
endif
ifneq ($(NOBLKMEM),$(LAST_NOBLKMEM))
		@-touch -c $(ALLSRC)
endif
ifneq ($(NOBUFMEM),$(LAST_NOBUFMEM))
		@-touch -c $(ALLSRC)
endif
ifneq ($(NOBLKBUFMEM),$(LAST_NOBLKBUFMEM))
		@-touch -c $(ALLSRC)
endif
		@-rm -f $(LASTSETTINGS)
		@echo "LAST_PARASCIP=$(PARASCIP)" >> $(LASTSETTINGS)
		@echo "LAST_SHARED=$(SHARED)" >> $(LASTSETTINGS)
		@echo "LAST_USRFLAGS=$(USRFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USROFLAGS=$(USROFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRCFLAGS=$(USRCFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRCXXFLAGS=$(USRCXXFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRLDFLAGS=$(USRLDFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRARFLAGS=$(USRARFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRDFLAGS=$(USRDFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_NOBLKMEM=$(NOBLKMEM)" >> $(LASTSETTINGS)
		@echo "LAST_NOBUFMEM=$(NOBUFMEM)" >> $(LASTSETTINGS)
		@echo "LAST_NOBLKBUFMEM=$(NOBLKBUFMEM)" >> $(LASTSETTINGS)

$(LINKSMARKERFILE):
		@$(MAKE) links

.PHONY: links
links:		| $(SCIPSDPLIBDIR) echosoftlinks $(SOFTLINKS)
		@rm -f $(LINKSMARKERFILE)
		@echo "this is only a marker" > $(LINKSMARKERFILE)

.PHONY: echosoftlinks
echosoftlinks:
		@echo
		@echo "- Current settings: SDPS=$(SDPS) LPS=$(LPS) SUFFIX=$(LINKLIBSUFFIX) OSTYPE=$(OSTYPE) ARCH=$(ARCH) COMP=$(COMP)"
		@echo
		@echo "* SCIPSDP needs some softlinks to external programs, in particular, SDP-solvers."
		@echo "* Please insert the paths to the corresponding directories/libraries below."
		@echo "* The links will be installed in the 'lib' directory."
		@echo "* For more information and if you experience problems see the INSTALL file."
		@echo
		@echo -e $(SDPIINSTMSG)

.PHONY: $(SOFTLINKS)
$(SOFTLINKS):
ifeq ($(MAKESOFTLINKS), true)
		@$(SHELL) -ec 'if test ! -e $@ ; \
			then \
				DIRNAME=`dirname $@` ; \
				BASENAMEA=`basename $@ .$(STATICLIBEXT)` ; \
				BASENAMESO=`basename $@ .$(SHAREDLIBEXT)` ; \
				echo ; \
				echo "- preparing missing soft-link \"$@\":" ; \
				if test -e $$DIRNAME/$$BASENAMEA.$(SHAREDLIBEXT) ; \
				then \
					echo "* this soft-link is not necessarily needed since \"$$DIRNAME/$$BASENAMEA.$(SHAREDLIBEXT)\" already exists - press return to skip" ; \
				fi ; \
				if test -e $$DIRNAME/$$BASENAMESO.$(STATICLIBEXT) ; \
				then \
					echo "* this soft-link is not necessarily needed since \"$$DIRNAME/$$BASENAMESO.$(STATICLIBEXT)\" already exists - press return to skip" ; \
				fi ; \
				echo "> Enter soft-link target file or directory for \"$@\" (return if not needed): " ; \
				echo -n "> " ; \
				cd $$DIRNAME ; \
				eval $(READ) TARGET ; \
				cd $(SCIPSDPDIR) ; \
				if test "$$TARGET" != "" ; \
				then \
					echo "-> creating softlink \"$@\" -> \"$$TARGET\"" ; \
					rm -f $@ ; \
					$(LN_s) $$TARGET $@ ; \
				else \
					echo "* skipped creation of softlink \"$@\". Call \"make links\" if needed later." ; \
				fi ; \
				echo ; \
			fi'
endif

#-----------------------------------------------------------------------------
.PHONY: test
test:
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) \
		$(DISPFREQ) $(CONTINUE) $(LOCK) $(SCIPSDPVERSION) $(SDPS) $(VALGRIND) $(CLIENTTMPDIR) $(OPTCOMMAND);

# include local targets
-include make/local/make.targets

.PHONY: testcluster
testcluster:
		cd check; \
		$(SHELL) ./check_cluster.sh $(TEST) $(MAINFILE) $(SETTINGS) \
		$(notdir $(MAINFILE)) $(TIME) $(NODES) $(MEM) \
		$(THREADS) $(FEASTOL) $(LPS) $(DISPFREQ) $(CONTINUE) \
		$(QUEUETYPE) $(QUEUE) $(PPN) $(CLIENTTMPDIR) \
		$(NOWAITCLUSTER) $(EXCLUSIVE) $(PERMUTE) $(OPTCOMMAND);

#-----------------------------------------------------------------------------

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(SDPIINC) $(DFLAGS) $(MAINCCSRC) $(SDPICCSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/scipsdp/\([0-9A-Za-z\_]*\).cpp|$$\(OBJDIR\)/\2.o: $(SRCDIR)/scipsdp/\2.cpp|g'\'' \
		>$(MAINDEP)'
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(SDPIINC) $(DFLAGS) $(MAINCSRC) $(SDPICSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/scipsdp/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/scipsdp/\2.c|g'\'' \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/sdpi/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/sdpi/\2.c|g'\'' \
		>>$(MAINDEP)'

-include	$(MAINDEP)

$(MAINFILE):	preprocess $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINCOBJFILES) $(MAINCCOBJFILES) $(SDPIOBJ) | $(SDPOBJSUBDIRS) $(BINDIR)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINCOBJFILES) $(MAINCCOBJFILES) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
                $(OFLAGS) $(LPSLDFLAGS) \
		$(SDPIOBJ) $(SDPILIB) $(LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c | $(SDPOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(SDPIINC) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp | $(SDPOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(SDPIINC) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
