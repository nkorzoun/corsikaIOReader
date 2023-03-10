#=============================================================================
#    corsikaIOreader is a tool to read CORSIKA eventio files
#    Copyright (C) 2004, 2013, 2019 Gernot Maier and Henrike Fleischhack
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#=============================================================================
#
# Makefile for corsikaIOreader
#

ARCH := $(shell uname)

# linux flags
ifeq ($(ARCH),Linux)
CXX           = g++ 
CXXFLAGS      = -g -O3 -Wall -fPIC -fno-strict-aliasing  -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -D_LARGEFILE64_SOURCE
LD            = g++
LDFLAGS       = -O
# LDFLAGS       =  -pg -O
SOFLAGS       = -shared
endif

# Apple OS X flags
ifeq ($(ARCH),Darwin)
CXX           = clang++ 
CXXFLAGS      = -g -O3 -Wall -fPIC  -fno-strict-aliasing
LD            = clang++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

OutPutOpt     = -o

# Root

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)
LIBS         += -lMinuit
LIBS 		 += -lTreePlayer -lGenVector -lGeomBuilder -lGeomPainter -lGeom -lGed  #NK

GLIBS         = $(ROOTGLIBS)
GLIBS        += -lMinuit
GLIBS 		 += -lTreePlayer -lGenVector -lGeomBuilder -lGeomPainter -lGeom -lGed  #NK

INCLUDEFLAGS  = -I. -I./inc/
CXXFLAGS     += $(INCLUDEFLAGS)

vpath %.h ./src/ ./inc
vpath %.cpp ./src/
vpath %.c ./src/

.cpp.o:
	$(CXX) $(CXXFLAGS)  -c $<

all:	corsikaIOreader


corsikaIOreader:	GUtilityFuncts.o GDefinition.o straux.o eventio.o warning.o io_simtel.o VIOHistograms.o atmo.o fileopen.o sim_cors.o corsikaIOreader.o VAtmosAbsorption.o VGrisu.o VCORSIKARunheader.o VCORSIKARunheader_Dict.o
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt) $@
		@echo "$@ done"

clean:	
	rm -f *.o *_Dict*

.SUFFIXES: .o

atmo.o: atmo.h fileopen.h
eventio.o: initial.h io_basic.h 
fileopen.o: initial.h straux.h fileopen.h
iact.o: initial.h io_basic.h mc_tel.h
io_simtel.o: initial.h io_basic.h mc_tel.h
corsikaIOreader.o: initial.h io_basic.h mc_tel.h atmo.h sim_cors.h
warning.o:	warning.h initial.h
straux.o: initial.h straux.h
VIOHistograms.o:	mc_tel.h sim_cors.h
VAtmosAbsorption.o:	VAtmosAbsorption.h
VCORSIKARunheader.o:	VCORSIKARunheader.h
VGrisu.o:	mc_tel.h sim_cors.h VCORSIKARunheader.h
sim_cors.o:	sim_cors.h
GUtilityFuncts.o: GUtilityFuncts.h #NK
GDefinition.o: GDefinition.h #NK

VCORSIKARunheader_Dict.cpp:	VCORSIKARunheader.h VCORSIKARunheaderLinkDef.h
	@echo "Generating dictionary $@..."
	@echo rootcint -f $@  -c -p inc/VCORSIKARunheader.h $^
	@rootcint -f $@  -c -p inc/VCORSIKARunheader.h $^ 

