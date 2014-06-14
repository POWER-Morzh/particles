# This is a makefile for Peano's PoissonWithJacobi component
# generated by the Peano Development Tools (PDT) 

# Include files
-include files.mk

MANTISSA ?= 8
MUL_FACTOR ?= 10



# Set Paths
# ---------
# Please adopt these directories to your Peano installation. The values are 
# used by the compiler and linker calls.
PEANO_HOME   = ../../src

# this is the path/root where your projects are stored
PROJECT_ROOT = ..

# points to directory holding the project's sources
PROJECT_HOME = particles  


# Set Dimension
# -------------
DIM=-DDim2
#DIM=-DDim3
#DIM=-DDim4


# Configure Peano
#----------------
PROJECT_CFLAGS = #-DDebug -DAsserts 
PROJECT_LFLAGS = 
EXECUTABLE=bin/peano-Particles-F$(MUL_FACTOR)-M$(MANTISSA)
#EXECUTABLE=peano-PoissonWithJacobi-debug


# Configure System
# ----------------
# These values are used to assemble the symbols SYSTEM_CFLAGS and SYSTEM_LFLAGS.
INCLUDE_TBB=$(TBB_INC)
#INCLUDE_OMP=$(OMP_INC)
INCLUDE_OMP=-fopenmp
INCLUDE_MPI=$(MPI_INC)


LINK_TBB=$(TBB_SHLIB)
#LINK_OMP=$(OMP_SHLIB)
LINK_OMP=-fopenmp
#LINK_MPI=-lpthread -lrt -lmpich
LINK_MPI=


# Assemble Compiler Flags
# -----------------------
SYSTEM_CFLAGS = $(INCLUDE_TBB) $(INCLUDE_MPI)
#SYSTEM_CFLAGS =  $(INCLUDE_OMP) $(INCLUDE_MPI)
SYSTEM_LFLAGS = $(LINK_TBB)    $(LINK_MPI)
#SYSTEM_LFLAGS =  $(LINK_OMP)    $(LINK_MPI)



# Settings for the GNU Compiler (Debug Mode)
# ------------------------------------------
CC=g++
COMPILER_CFLAGS=-O0  -Wstrict-aliasing -fstrict-aliasing -ggdb -std=c++0x -DMANTISSA=$(MANTISSA) -DMUL_FACTOR=$(MUL_FACTOR)
#-pedantic -pedantic-errors -Wall -Werror

# Settings for the GNU Compiler (Release Mode)
# --------------------------------------------
#CC=g++
#COMPILER_CFLAGS=-O3 -fstrict-aliasing -fno-rtti -fno-exceptions -std=c++0x


# Settings for the Intel Compiler (Debug Mode)
# --------------------------------------------
#CC=icpc
#COMPILER_CFLAGS=-O0 -fstrict-aliasing -std=c++0x


# Settings for the Intel Compiler (Release Mode)
# ----------------------------------------------
#CC=icpc
#COMPILER_CFLAGS=-fast -fstrict-aliasing  -DMPICH_IGNORE_CXX_SEEK -fast -fno-rtti -no-ipo -ip -xSSE4.2 -std=c++0x


# Settings for the XLC Compiler (Debug Mode)
# --------------------------------------------
#CC=xlcc
#COMPILER_CFLAGS=-O0 -fstrict-aliasing -qpack_semantic=gnu 


# Settings for the XLC Compiler (Release Mode)
# ----------------------------------------------
#CC=xlcc
#COMPILER_CFLAGS=-fast -fstrict-aliasing -qpack_semantic=gnu 



#
# Linker Settings
# ---------------
# By default, I use the compiler command. But you might wanna modify it.
LL=$(CC)
COMPILER_LFLAGS=


OBJECTS=$(SOURCES:.cpp=.o)



all: header build
	@echo "value is $(MANTISSA)"


files.mk:
	touch files.mk
	echo -n SOURCES= > files.mk
	find $(PEANO_HOME) -name '*.cpp' -printf '%p ' >> files.mk
	find $(PROJECT_ROOT)/$(PROJECT_HOME) -name '*.cpp' -printf '%p ' >> files.mk
	find $(PEANO_HOME) -name '*.cpp' -printf 'COMPILE %p TO %p OBJECT FILE\n' >> compiler-minutes.txt
	find $(PROJECT_ROOT)/$(PROJECT_HOME) -name '*.cpp' -printf 'COMPILE %p TO %p OBJECT FILE\n' >> compiler-minutes.txt
	echo -n '\n\nLINK OBJECT FILES OF ' >> compiler-minutes.txt
	find $(PEANO_HOME) -name '*.cpp' -printf '%p ' >> compiler-minutes.txt
	find $(PROJECT_ROOT)/$(PROJECT_HOME) -name '*.cpp' -printf '%p ' >> compiler-minutes.txt



header:
	@echo  --- This is Peano 3 ---


build: files.mk $(OBJECTS)
	$(LL) $(OBJECTS) -o $(EXECUTABLE) $(PROJECT_LFLAGS) $(COMPILER_LFLAGS) $(SYSTEM_LFLAGS) 
	@echo
	@echo build of Peano with component Particles successful
	@echo visit Peano\'s homepage at http://www.peano-framework.org


clean:
	#rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)
	rm -f files.mk
	rm -f compiler-minutes.txt


$(OBJECTS): %.o : %.cpp
	$(CC) $(DIM) $(PROJECT_CFLAGS) $(COMPILER_CFLAGS) $(SYSTEM_CFLAGS) -I$(PROJECT_ROOT) -I$(PEANO_HOME) -c $< -o $@
