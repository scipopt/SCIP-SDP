#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*                                                                           */
#/*                  This file is part of the                                 */
#/*      SDP-Package for SCIP: a solving framework for                        */
#/*                            mixed-integer semidefinite programms           */
#/*                                                                           */
#/* Copyright (C) 2011-2014 Discrete Optimization, TU Darmstadt               */
#/*                         EDOM, FAU Erlangen-Nürnberg                       */
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
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#@file    Makefile
#@brief   Makefile for C++ SDP-Interface for SCIP
#@author  Sonja Mars, Lars Schewe, Marc Pfetsch, Tristan Gally

#-----------------------------------------------------------------------------
# own variables
#-----------------------------------------------------------------------------

SCIPSDPVERSION	=	1.0
SDPS		=	none

GCCWARN		+= 	-Wextra

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
FLAGS 		+= 	-I$(DSDP_INCLUDE_DIR) -DUSE_DSDP
endif


SDPOBJSUBDIRS	=	$(OBJDIR)/scipsdp \
			$(OBJDIR)/sdpi


#-----------------------------------------------------------------------------
# main program
#-----------------------------------------------------------------------------

MAINNAME=	scipsdp
MAINOBJ	=	scipsdp/main.o \
		scipsdp/relax_sdp.o \
		scipsdp/objconshdlr_sdp.o \
		scipsdp/objreader_sdpa.o \
		scipsdp/SdpVarMapper.o \
		scipsdp/SdpCone.o \
		scipsdp/ScipStreamBuffer.o \
		scipsdp/SdpProblem.o 

ifeq ($(SDPS), dsdp)
MAINOBJ 	+= 	sdpi/sdpi_dsdp.o
endif

MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINDEP		=	$(SRCDIR)/depend.cppmain.$(OPT)

# @todo possibly add LPS
MAINFILE	=	$(BINDIR)/$(MAINNAME).$(BASE).$(SDPS)$(EXEEXTENSION)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))


#-----------------------------------------------------------------------------
# rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: tags
tags:
		rm -f TAGS; ctags -e src/*/*.c src/*/*.cpp src/*/*.h $(SCIPDIR)/src/*/*.c $(SCIPDIR)/src/*/*.h;

.PHONY: lint
lint:		$(MAINSRC)
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
		-rm -f $(OBJDIR)/*.o
		-rmdir $(OBJDIR)
endif
		-rm -f $(MAINFILE)

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).cpp|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.cpp|g'\'' \
		>$(MAINDEP)'

-include	$(MAINDEP)

$(MAINFILE):	$(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES) | $(SDPOBJSUBDIRS) $(BINDIR)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) \
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
