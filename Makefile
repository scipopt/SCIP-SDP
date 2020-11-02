#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*                                                                           */
#/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
#/* semidefinite programs based on SCIP.                                      */
#/*                                                                           */
#/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
#/*                         EDOM, FAU Erlangen-Nürnberg                       */
#/*               2014-2018 Discrete Optimization, TU Darmstadt               */
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
#/* Copyright (C) 2002-2018 Zuse Institute Berlin                             */
#/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
#/* see file COPYING in the SCIP distribution.                                */
#/*                                                                           */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#@file    Makefile
#@brief   Makefile for C++ SDP-Interface for SCIP


SCIPSDPDIR	=	./
# mark that this is a SCIPSDP internal makefile
SCIPSDPINTERNAL	=	true

# load project makefile
include $(SCIPSDPDIR)/make/make.scipsdpproj

SCIPREALPATH	=	$(realpath $(SCIPDIR))

# overwrite flags for dependencies
DFLAGS		=       -MM

#-----------------------------------------------------------------------------
# DSDP solver
SDPIOPTIONS	=	dsdp
ifeq ($(SDPS),dsdp)
SDPIINC		= 	-I$(SCIPSDPLIBDIR)/include/dsdpinc
SDPICSRC 	= 	src/sdpi/sdpisolver_dsdp.c
SDPIOBJ 	= 	$(OBJDIR)/sdpi/sdpisolver_dsdp.o
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/include/dsdpinc
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libdsdp.$(STATICLIBEXT)
SDPIINSTMSG	=	" -> \"dsdpinc\" is the path to the DSDP \"include\" directory, e.g., \"<DSDP-path>/include\".\n"
SDPIINSTMSG	+=	" -> \"libdsdp.*\" is the path to the DSDP library, e.g., \"<DSDP-path>/lib/libdsdp.$(STATICLIBEXT)\""
endif

#-----------------------------------------------------------------------------
# SDPA solver version >= 7.3.8
SDPIOPTIONS	+=	sdpa
ifeq ($(SDPS),sdpa)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/include/sdpainc
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libsdpa.$(STATICLIBEXT)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/shared/libsdpa.$(SHAREDLIBEXT)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/include/mumpsinc
ifeq ($(MUMPSSEQ),true)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libdmumps_seq.$(STATICLIBEXT)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libmumps_common_seq.$(STATICLIBEXT)
else
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libdmumps.$(STATICLIBEXT)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libmumps_common.$(STATICLIBEXT)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libpord.$(STATICLIBEXT)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libmpiseq.$(STATICLIBEXT)
endif
ifeq ($(OPENBLAS),true)
ifeq ($(SHARED),true)
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/shared/libopenblas.$(SHAREDLIBEXT).0
else
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/static/libopenblas.$(STATICLIBEXT)
endif
endif
SDPIINSTMSG	=	" -> \"sdpainc\" is the path to the SDPA \"include\" directory, e.g., \"<SDPA-path>/include\".\n"
SDPIINSTMSG	+=	" -> \"libsdpa.*\" is the path to the SDPA library, e.g., \"<SDPA-path>/lib/libsdpa.a\".\n"
SDPIINSTMSG	+=	" -> \"mumpsinc\" is the path to the mumps \"include\" directory, e.g., \"<SDPA-path>/mumps/include\".\n"
SDPIINSTMSG	+=	" -> \"libdmumps.*\" is the path to the dmumps library, e.g., \"<SDPA-path>/mumps/build/lib/libdmumps.$(STATICLIBEXT)\".\n"
SDPIINSTMSG	+=	" -> \"libdmumps_common.*\" is the path to the mumps_common library, e.g., \"<SDPA-path>/mumps/build/lib/libmumps_common.$(STATICLIBEXT)\".\n"
SDPIINSTMSG	+=	" -> \"libpord.*\" is the path to the pord library, e.g., \"<SDPA-path>/mumps/build/lib/libpord.$(STATICLIBEXT)\".\n"
SDPIINSTMSG	+=	" -> \"libmpiseq.*\" is the path to the mpiseq library, e.g., \"<SDPA-path>/mumps/build/libseq/libmpiseq.$(STATICLIBEXT)\".\n"
ifeq ($(OPENBLAS),true)
ifeq ($(SHARED),true)
SDPIINSTMSG	+=	" -> \"libopenblas.$(SHAREDLIBEXT).0\" is the openblas library.\n"
else
SDPIINSTMSG	+=	" -> \"libopenblas.$(STATICLIBEXT)\" is the openblas library.\n"
endif
endif
SDPIINC		=  	-I$(SCIPSDPLIBDIR)/include/sdpainc
SDPIINC		+= 	-I$(SCIPSDPLIBDIR)/include/mumpsinc
SDPICCSRC 	= 	src/sdpi/sdpisolver_sdpa.cpp
SDPICSRC	=
SDPIOBJ 	= 	$(OBJDIR)/sdpi/sdpisolver_sdpa.o

# disable OMP
ifeq ($(MUMPSSEQ),true)
OMP 		= 	false
endif
endif

#-----------------------------------------------------------------------------
# MOSEK solver
SDPIOPTIONS	+=	msk
ifeq ($(SDPS),msk)
# decide on 32 or 64 bit
BITEXT     	=  	$(word 2, $(subst _, ,$(ARCH)))
SDPIINC		= 	-I$(SCIPSDPLIBDIR)/include/mosekh
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/include/mosekh
SOFTLINKS	+=	$(SCIPSDPLIBDIR)/shared/libmosek$(BITEXT).$(SHAREDLIBEXT)
SDPIINSTMSG	=	"  -> \"mosekh\" is the path to the MOSEK \"h\" directory, e.g., \"<MOSEK-path>/8/tools/platform/linux64x86/h\".\n"
SDPIINSTMSG	+=	" -> \"libmosek$(BITEXT).*\" is the path to the MOSEK library, e.g., \"<MOSEK-path>/8/tools/platform/linux64x86/bin/libmosek$(BITEXT).$(SHAREDLIBEXT)\".\n"
SDPICSRC 	= 	src/sdpi/sdpisolver_mosek.c
SDPIOBJ 	= 	$(OBJDIR)/sdpi/sdpisolver_mosek.o
endif

#-----------------------------------------------------------------------------
# no solver
SDPIOPTIONS	+=	none
ifeq ($(SDPS),none)
SDPILIB		= 	-L$(SCIPSDPLIBDIR) -llapack -lblas
SDPICSRC 	= 	src/sdpi/sdpisolver_none.c
SDPIOBJ 	= 	$(OBJDIR)/sdpi/sdpisolver_none.o
SETTINGS	= 	lp_approx
endif

LINKSMARKERFILE	=	$(LIBDIR)/linkscreated.$(LPS)-$(LPSOPT).$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX).$(ZIMPL)-$(ZIMPLOPT).$(IPOPT)-$(IPOPTOPT).$(GAMS)


#-----------------------------------------------------------------------------

SDPOBJSUBDIRS	=	$(OBJDIR)/scipsdp \
			$(OBJDIR)/sdpi

#-----------------------------------------------------------------------------
# OMPSETTINGS (used to set number of threads for Openblas)
#-----------------------------------------------------------------------------

ifeq ($(OMP),true)
OMPFLAGS += -DOMP
endif

# handle long long int option
ifeq ($(LAPACKLONGINT),true)
CFLAGS		+=	-DLAPACK_LONGLONGINT
endif

#-----------------------------------------------------------------------------
# SCIPSDP
#-----------------------------------------------------------------------------

SCIPSDPNAME	=	scipsdp
SCIPSDPCOBJ	=	scipsdp/SdpVarmapper.o \
			scipsdp/SdpVarfixer.o \
			scipsdp/cons_sdp.o \
			scipsdp/cons_savedsdpsettings.o \
			scipsdp/cons_savesdpsol.o \
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
			scipsdp/reader_cbf.o \
			scipsdp/reader_sdpa_firsttry.o \
			scipsdp/prop_sdpobbt.o \
			scipsdp/prop_companalcent.o \
			scipsdp/table_relaxsdp.o \
			scipsdp/table_slater.o \
			scipsdp/table_sdpsolversuccess.o \
			sdpi/sdpi.o \
			sdpi/sdpsolchecker.o \
			sdpi/lapack_interface.o \
			scipsdpgithash.o

SCIPSDPCCOBJ 	=	scipsdp/scipsdpdefplugins.o \
			scipsdp/ScipStreamBuffer.o

SCIPSDPCSRC	=	$(addprefix $(SRCDIR)/,$(SCIPSDPCOBJ:.o=.c))
SCIPSDPCCSRC 	=	$(addprefix $(SRCDIR)/,$(SCIPSDPCCOBJ:.o=.cpp))
SCIPSDPDEP 	=	$(SRCDIR)/depend.cppmain.$(OPT)

SCIPSDPGITHASHFILE	= 	$(SRCDIR)/scipsdpgithash.c

# @todo possibly add LPS
SCIPSDPBIN		=	$(BINDIR)/$(SCIPSDPNAME).$(BASE).$(SDPS)$(EXEEXTENSION)
SCIPSDPSHORTLINK	=	$(BINDIR)/$(SCIPSDPNAME)
SCIPSDPCOBJFILES	=	$(addprefix $(OBJDIR)/,$(SCIPSDPCOBJ))
SCIPSDPCCOBJFILES	=	$(addprefix $(OBJDIR)/,$(SCIPSDPCCOBJ))

MAINOBJ		=	scipsdp/main.o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.c))
MAINOBJFILES  	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

ALLSRC		=	$(SCIPSDPCSRC) $(SCIPSDPCCSRC) $(SDPICSRC) $(SDPICCSRC) $(MAINSRC)
LINKSMARKERFILE =	$(SCIPSDPLIBDIR)/linkscreated.$(SDPS).$(LPS)-$(LPSOPT).$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
LASTSETTINGS 	=	$(OBJDIR)/make.lastsettings

SCIPSDPLIBOBJFILES	=	$(addprefix $(OBJDIR)/,$(SCIPSDPCOBJ))
SCIPSDPLIBOBJFILES	+=	$(addprefix $(OBJDIR)/,$(SCIPSDPCCOBJ))
SCIPSDPLIBOBJFILES	+=	$(SDPIOBJ)

#-----------------------------------------------------------------------------
# rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(SCIPSDPBIN) $(SCIPSDPCOBJFILES) $(SCIPSDPCCOBJFILES) $(MAINOBJFILES) $(SCIPSDPLIBFILE) $(SDPIOBJ) $(SCIPSDPSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(SCIPSDPBIN) $(SCIPSDPSHORTLINK)

.PHONY: checkdefines
checkdefines:
ifeq ($(SDPIOBJ),)
		$(error invalid SDP solver selected: SDPS=$(SDPS). Possible options are: $(SDPIOPTIONS))
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
		rm -f TAGS; ctags -e src/*/*.c src/*/*.cpp src/*/*.h $(SCIPREALPATH)/src/*/*.c $(SCIPREALPATH)/src/*/*.h; sed -i 's!\#undef .*!!g' TAGS

# include target to detect the current git hash
-include make/make.detectgithash

# this empty target is needed for the SCIP release versions
githash::      # do not remove the double-colon

.PHONY: lint
lint:		$(SCIPSDPCSRC) $(SCIPSDPCCSRC) $(SDPICSRC) $(SDPICCSRC) $(MAINSRC)
		-rm -f lint.out
ifeq ($(FILES),)
	$(SHELL) -ec 'for i in $^; \
			do \
				echo $$i; \
				$(LINT) -I$(SCIPREALPATH) lint/main-gcc.lnt +os\(lint.out\) -u -zero \
				$(USRFLAGS) $(FLAGS) $(SDPIINC) -I/usr/include -UNDEBUG -USCIP_WITH_READLINE -USCIP_ROUNDING_FE -D_BSD_SOURCE $$i; \
			done'
else
		$(SHELL) -ec  'for i in $(FILES); \
			do \
				echo $$i; \
				$(LINT) -I$(SCIPREALPATH) lint/main-gcc.lnt +os\(lint.out\) -u -zero \
				$(USRFLAGS) $(FLAGS) $(SDPIINC) -I/usr/include -UNDEBUG -USCIP_WITH_READLINE -USCIP_ROUNDING_FE -D_BSD_SOURCE $$i; \
			done'
endif

.PHONY: pclint
pclint:		$(SCIPSDPCSRC) $(SCIPSDPCCSRC) $(SDPICSRC) $(SDPICCSRC) $(MAINSRC)
		-rm -f pclint.out
ifeq ($(FILES),)
		@$(SHELL) -ec 'echo "-> running pclint ..."; \
			for i in $^; \
			do \
				echo $$i; \
				$(PCLINT) -I$(SCIPREALPATH) -I$(SCIPREALPATH)/pclint main-gcc.lnt +os\(pclint.out\) -b -u -zero \
				$(USRFLAGS) $(FLAGS) $(SDPIINC) -uNDEBUG -uSCIP_WITH_READLINE -uSCIP_ROUNDING_FE -D_BSD_SOURCE $$i; \
			done'
else
		@$(SHELL) -ec  'echo "-> running pclint on specified files ..."; \
			for i in $(FILES); \
			do \
				echo $$i; \
				$(PCLINT) -I$(SCIPREALPATH) -I$(SCIPREALPATH)/pclint main-gcc.lnt +os\(pclint.out\) -b -u -zero \
				$(USRFLAGS) $(FLAGS) $(SDPIINC) -uNDEBUG -uSCIP_WITH_READLINE -uSCIP_ROUNDING_FE -D_BSD_SOURCE $$i; \
			done'
endif


.PHONY: doc
doc:
		cd doc; $(DOXY) $(SCIPSDPNAME).dxy

$(SCIPSDPSHORTLINK): $(SCIPSDPBIN)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(SCIPSDPBIN)) $(notdir $@)

$(OBJDIR):
		@mkdir -p $(OBJDIR);

$(SDPOBJSUBDIRS):	| $(OBJDIR)
		@-mkdir -p $(SDPOBJSUBDIRS);

$(SCIPSDPLIBDIR):
		@-mkdir -p $(SCIPSDPLIBDIR)

$(SCIPSDPLIBDIR)/include: $(SCIPSDPLIBDIR)
		@-mkdir -p $(SCIPSDPLIBDIR)/include

$(SCIPSDPLIBDIR)/static: $(SCIPSDPLIBDIR)
		@-mkdir -p $(SCIPSDPLIBDIR)/static

$(SCIPSDPLIBDIR)/shared: $(SCIPSDPLIBDIR)
		@-mkdir -p $(SCIPSDPLIBDIR)/shared

$(BINDIR):
		-@test -d $(BINDIR) || { \
		echo "-> Creating $(BINDIR) directory"; \
		mkdir -p $(BINDIR); }

.PHONY: libscipsdp
libscipsdp:	preprocess
		@$(MAKE) $(SCIPSDPLIBFILE) $(SCIPSDPLIBLINK) $(SCIPSDPLIBSHORTLINK)

$(SCIPSDPLIBFILE):	$(SCIPSDPLIBOBJFILES) | $(SCIPSDPLIBDIR)/$(LIBTYPE)
		@echo "-> generating library $@"
		-rm -f $@
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(SCIPSDPLIBOBJFILES)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

.PHONY: clean
clean:
ifneq ($(OBJDIR),)
		@-rm -f $(LASTSETTINGS)
		@-rm -f $(OBJDIR)/scipsdp/*.o
		@-rm -f $(OBJDIR)/sdpi/*.o
		@-rm -f $(OBJDIR)/*.o
		@-rmdir $(OBJDIR)/scipsdp
	 	@-rmdir $(OBJDIR)/sdpi
		@-rmdir $(OBJDIR)
endif
		-rm -f $(SCIPSDPBIN)

#-----------------------------------------------------------------------------
-include $(LASTSETTINGS)

.PHONY: touchexternal
touchexternal:	$(SCIPSDPLIBDIR) $(OBJDIR)
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
ifneq ($(SDPS),$(LAST_SDPS))
		@-touch -c $(SRCDIR)/scipsdp/cons_sdp.c
endif
ifneq ($(OMP),$(LAST_OMP))
		@-touch -c $(SRCDIR)/scipsdp/cons_sdp.c
endif
ifneq ($(SCIPSDPGITHASH),$(LAST_SCIPSDPGITHASH))
		@$(MAKE) githash
endif
		@$(SHELL) -ec 'if test ! -e $(SCIPSDPGITHASHFILE) ; \
			then \
				echo "-> generating $(SCIPSDPGITHASHFILE)" ; \
				$(MAKE) githash ; \
			fi'
		@-rm -f $(LASTSETTINGS)
		@echo "LAST_PARASCIP=$(PARASCIP)" >> $(LASTSETTINGS)
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
		@echo "LAST_SDPS=$(SDPS)" >> $(LASTSETTINGS)
		@echo "LAST_OMP=$(OMP)" >> $(LASTSETTINGS)
		@echo "LAST_SCIPSDPGITHASH=$(SCIPSDPGITHASH)" >> $(LASTSETTINGS)

$(LINKSMARKERFILE):
		@$(MAKE) links

.PHONY: links
links:		| $(SCIPSDPLIBDIR) $(SCIPSDPLIBDIR)/include $(SCIPSDPLIBDIR)/static $(SCIPSDPLIBDIR)/shared echosoftlinks $(SOFTLINKS)
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
				cd $(SCIPSDPREALPATH) ; \
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
		@-(cd check && ln -fs $(SCIPREALPATH)/check/check.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/allcmpres.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/cmpres.awk);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/evalcheck.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/evalcheck_cluster.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/evaluate.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/check.awk);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/cmpres.awk);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/getlastprob.awk);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/configuration_solufile.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/configuration_set.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/configuration_logfiles.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/configuration_tmpfile_setup_scip.sh configuration_tmpfile_setup_$(SCIPSDPNAME).sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/run.sh);
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(SCIPSDPBIN) $(SETTINGS) $(notdir $(SCIPSDPBIN)) $(OUTPUTDIR) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) \
			$(CONTINUE) $(LOCK) $(SCIPSDPVERSION) $(SDPS) $(DEBUGTOOL) $(CLIENTTMPDIR) $(REOPT) $(OPTCOMMAND) $(SETCUTOFF) $(MAXJOBS) $(VISUALIZE) \
			$(PERMUTE) $(SEEDS) $(GLBSEEDSHIFT) $(STARTPERM);

# include local targets
-include make/local/make.targets

.PHONY: testcluster
testcluster:
		@-(cd check && ln -fs $(SCIPREALPATH)/check/check.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/allcmpres.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/cmpres.awk);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/evalcheck.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/evalcheck_cluster.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/evaluate.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/check_cluster.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/check.awk);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/cmpres.awk);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/getlastprob.awk);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/configuration_solufile.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/configuration_cluster.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/configuration_set.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/configuration_logfiles.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/configuration_tmpfile_setup_scip.sh configuration_tmpfile_setup_$(SCIPSDPNAME).sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/run.sh);
		@-(cd check && ln -fs $(SCIPREALPATH)/check/testfiles.sh);
		cd check; \
		$(SHELL) ./check_cluster.sh $(TEST) $(PWD)/$(SCIPSDPBIN) $(SETTINGS) $(notdir $(SCIPSDPBIN)) $(OUTPUTDIR) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(SDPS) $(DISPFREQ) \
			$(CONTINUE) $(QUEUETYPE) $(QUEUE) $(PPN) $(CLIENTTMPDIR) $(NOWAITCLUSTER) $(EXCLUSIVE) $(PERMUTE) $(SEEDS) $(GLBSEEDSHIFT) $(STARTPERM) $(DEBUGTOOL) $(REOPT) $(OPTCOMMAND) \
			$(SETCUTOFF) $(VISUALIZE) $(CLUSTERNODES) $(SLURMACCOUNT);

#-----------------------------------------------------------------------------

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCXX) $(DFLAGS) $(FLAGS) $(SDPIINC) $(DFLAGS) $(SCIPSDPCCSRC) $(SDPICCSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/scipsdp/\([0-9A-Za-z\_]*\).cpp|$$\(OBJDIR\)/\2.o: $(SRCDIR)/scipsdp/\2.cpp|g'\'' \
		>$(SCIPSDPDEP)'
		$(SHELL) -ec '$(DCXX) $(DFLAGS) $(FLAGS) $(SDPIINC) $(DFLAGS) $(SCIPSDPCSRC) $(SDPICSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/scipsdp/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/scipsdp/\2.c|g'\'' \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/sdpi/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/sdpi/\2.c|g'\'' \
		>>$(SCIPSDPDEP)'

-include	$(SCIPSDPDEP)

$(SCIPSDPBIN):	$(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) libscipsdp $(MAINOBJFILES) | $(SDPOBJSUBDIRS) $(BINDIR)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) -L$(SCIPSDPLIBDIR)/$(LIBTYPE) -l$(SCIPSDPLIB) $(SDPILIB) $(LINKCXXSCIPALL) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c | $(SDPOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(SDPIINC) $(BINOFLAGS) $(CFLAGS) $(OMPFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp | $(SDPOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(SDPIINC) $(BINOFLAGS) $(CXXFLAGS) $(OMPFLAGS) -c $< $(CXX_o)$@


.PHONY: help
help:
		@echo "Use the SCIP-SDP makefile system."
		@echo
		@echo "  All options of the SCIP makefile apply here as well."
		@echo
		@echo "  Additional SCIP-SDP options:"
		@echo "  - SDPS={msk|dsdp|sdpa|sdpa740|none}: Determine SDP-solver."
		@echo "      msk: Mosek SDP-solver"
		@echo "      dsdp: DSDP SDP-solver"
		@echo "      sdpa: version >= 7.3.8 of the SDPA SDP-solver"
		@echo "      none: no SDP-solver"
		@echo "  - OPENBLAS={true|false}: use openblas"
		@echo "  - OMP={true|false}: use OMP"
		@echo "  - LAPACKLONGINT={true|false}: use long long ints for lapack (e.g., with SDPA 7.4.0 and openblas)"

#---- EOF --------------------------------------------------------------------
