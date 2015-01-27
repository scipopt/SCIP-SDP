#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*                                                                           */
#/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
#/* semidefinite programms based on SCIP.                                     */
#/*                                                                           */
#/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
#/*                         EDOM, FAU Erlangen-Nürnberg                       */
#/*               2014-2015 Discrete Optimization, TU Darmstadt               */
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
#/* Copyright (C) 2002-2015 Zuse Institute Berlin                             */
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

SCIPSDPVERSION	=	1.0
SDPS		=	none

GCCWARN		+= 	-Wextra

# check and testing parameters
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

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------

# possibly load local makefile
-include make.local

include $(SCIPDIR)/make/make.project


#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

LDFLAGS 	+= 	-lobjscip -llapack -lblas

ifeq ($(SDPS),dsdp)
LDFLAGS		+= 	-L$(DSDP_LIB_DIR) -ldsdp
FLAGS 		+= 	-I$(DSDP_INCLUDE_DIR)
endif

ifeq ($(SDPS),sdpa)
LDFLAGS         +=      -L$(SDPA_LIB_DIR) -lsdpa $(SDPA_LDFLAGS)
FLAGS           +=      -I$(SDPA_INCLUDE_DIR) $(SDPA_FLAGS)
endif

SDPOBJSUBDIRS	=	$(OBJDIR)/scipsdp \
			$(OBJDIR)/sdpi


#-----------------------------------------------------------------------------
# main program
#-----------------------------------------------------------------------------

MAINNAME=	scipsdp
MAINCOBJ=	scipsdp/SdpVarmapper.o \
		scipsdp/SdpVarfixer.o \
		scipsdp/disp_sdpiterations.o \
		scipsdp/disp_sdpavgiterations.o \
		scipsdp/prop_sdpredcost.o \
		scipsdp/branch_sdpmostfrac.o \
		scipsdp/branch_sdpmostinf.o \
		scipsdp/branch_sdpobjective.o \
		scipsdp/branch_sdpinfobjective.o \
	   scipsdp/branch_cs.o \
	   scipsdp/nodesel_prio.o \
		sdpi/sdpi.o

MAINCCOBJ=	scipsdp/main.o \
		scipsdp/relax_sdp.o \
		scipsdp/objreader_sdpa.o \
		scipsdp/cons_sdp.o \
		scipsdp/ScipStreamBuffer.o

ifeq ($(SDPS),dsdp)
MAINCOBJ 	+= 	sdpi/sdpisolver_dsdp.o
endif

ifeq ($(SDPS),sdpa)
MAINCOBJ 	+= 	sdpi/sdpisolver_sdpa.o
endif

MAINCSRC	=	$(addprefix $(SRCDIR)/,$(MAINCOBJ:.o=.c))
MAINCCSRC	=	$(addprefix $(SRCDIR)/,$(MAINCCOBJ:.o=.cpp))
MAINDEP		=	$(SRCDIR)/depend.cppmain.$(OPT)

# @todo possibly add LPS
MAINFILE	=	$(BINDIR)/$(MAINNAME).$(BASE).$(SDPS)$(EXEEXTENSION)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINCOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINCOBJ))
MAINCCOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINCCOBJ))


#-----------------------------------------------------------------------------
# rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINCOBJFILES) $(MAINCCOBJFILES) $(MAINSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: tags
tags:
		rm -f TAGS; ctags -e src/*/*.c src/*/*.cpp src/*/*.h $(SCIPDIR)/src/*/*.c $(SCIPDIR)/src/*/*.h;

.PHONY: lint
lint:		$(MAINCSRC) $(MAINCCSRC)
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

$(MAINSHORTLINK):	$(MAINFILE)
	@rm -f $@
	cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
	@mkdir -p $(OBJDIR);

$(SDPOBJSUBDIRS):	| $(OBJDIR)
	@-mkdir -p $(SDPOBJSUBDIRS);

$(BINDIR):
	-@test -d $(BINDIR) || { \
	echo "-> Creating $(BINDIR) directory"; \
	mkdir -p $(BINDIR); }

.PHONY: clean
clean:
ifneq ($(OBJDIR),)
		-rm -f $(OBJDIR)/scipsdp/*.o
		-rm -f $(OBJDIR)/sdpi/*.o
		-rmdir $(OBJDIR)/scipsdp
		-rmdir $(OBJDIR)/sdpi
		-rmdir $(OBJDIR)
endif
		-rm -f $(MAINFILE)

#-----------------------------------------------------------------------------
.PHONY: test
test:
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(SCIPSDPVERSION) $(SDPS) $(VALGRIND) $(CLIENTTMPDIR) $(OPTCOMMAND);

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
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(MAINCCSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/scipsdp/\([0-9A-Za-z\_]*\).cpp|$$\(OBJDIR\)/\2.o: $(SRCDIR)/scipsdp/\2.cpp|g'\'' \
		>$(MAINDEP)'
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(MAINCSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/scipsdp/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/scipsdp/\2.c|g'\'' \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/sdpi/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/sdpi/\2.c|g'\'' \
		>>$(MAINDEP)'

-include	$(MAINDEP)

$(MAINFILE):	$(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINCOBJFILES) $(MAINCCOBJFILES) | $(SDPOBJSUBDIRS) $(BINDIR)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINCOBJFILES) $(MAINCCOBJFILES) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
                $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c | $(SDPOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp | $(SDPOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
