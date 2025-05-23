#---------------------------------------------------------------------------
#
#  Makefile for building the Skin++ library and load module
#  Author : A. Hebert (2018-5-10)
#
#---------------------------------------------------------------------------
#
nbit = -m32
ARCH = $(shell uname -m)
ifeq (,$(filter $(ARCH),i386 i686))
  nbit = -m64
endif

DIRNAME = $(shell uname -sm | sed 's/[ ]/_/')
OS = $(shell uname -s | cut -d"_" -f1)
ccompiler = g++
opt = -O -g
ifeq ($(openmp),1)
  COMP = -fopenmp
  FOMP = -fopenmp -D_OPENMP
  ccompiler = g++-11
else
  COMP =
  FOMP =
endif

ifeq ($(intel),1)
  fcompiler = ifort
else
  fcompiler = gfortran
endif

clib = -lstdc++
ifeq ($(OS),Darwin)
  F90 = $(fcompiler)
  C = $(ccompiler)
  FLAGS = -DLinux -DUnix
  CFLAGS = -Wall $(nbit) -fPIC
  LFLAGS = $(nbit)
  clib = -lc++
else
ifeq ($(OS),Linux)
  F90 = $(fcompiler)
  C = $(ccompiler)
  FLAGS = -DLinux -DUnix
  CFLAGS = -Wall $(nbit) -fPIC
  LFLAGS = $(nbit)
else
ifeq ($(OS),CYGWIN)
  F90 = $(fcompiler)
  C = $(ccompiler)
  FLAGS = -DLinux -DUnix
  CFLAGS = -Wall $(nbit) -fPIC
  LFLAGS = $(nbit)
else
ifeq ($(OS),SunOS)
  MAKE = gmake
  F90 = f90
  C = cc
  FLAGS = -DSunOS -DUnix
  CFLAGS = $(nbit)
  LFLAGS = $(nbit)
else
ifeq ($(OS),AIX)
  opt = -O4
  MAKE = gmake
  DIRNAME = AIX
  F90 = xlf90
  C = xlc
  FLAGS = -DAIX -DUnix
  CFLAGS = -qstrict
  LFLAGS = -qstrict -bmaxdata:0x80000000 -qipa
else
  $(error $(OS) is not a valid OS)
endif
endif
endif
endif
endif
INCLUDE = -I../../Ganlib/src/

ifeq ($(intel),1)
  FFLAGS = -fPIC
  FFLAG77 = -fPIC
  lib = ../lib/$(DIRNAME)_intel
  bin = ../bin/$(DIRNAME)_intel
  libUtl = ../../Utilib/lib/$(DIRNAME)_intel
  libGan = ../../Ganlib/lib/$(DIRNAME)_intel
  libTri = ../../Trivac/lib/$(DIRNAME)_intel
  libDra = ../../Dragon/lib/$(DIRNAME)_intel
  libDon = ../../Donjon/lib/$(DIRNAME)_intel
else
  lib = ../lib/$(DIRNAME)
  bin = ../bin/$(DIRNAME)
  libUtl = ../../Utilib/lib/$(DIRNAME)
  libGan = ../../Ganlib/lib/$(DIRNAME)
  libTri = ../../Trivac/lib/$(DIRNAME)
  libDra = ../../Dragon/lib/$(DIRNAME)
  libDon = ../../Donjon/lib/$(DIRNAME)
endif

ifeq ($(hdf5),1)
  CFLAGS += -DHDF5_LIB
  LFLAGS += -L${HDF5_API} -lhdf5
endif

SRCC = $(shell ls *.cxx)
OBJC = $(SRCC:.cxx=.o)
all : sub-make Skin++
ifeq ($(openmp),1)
	@echo 'Skin++: openmp is defined'
endif
ifeq ($(intel),1)
	@echo 'Skin++: intel is defined'
endif
ifeq ($(hdf5),1)
	@echo 'Dragon: hdf5 is defined'
endif
sub-make:
	$(MAKE) openmp=$(openmp) intel=$(intel) hdf5=$(hdf5) -C ../../Donjon/src
%.o : %.cxx
	$(C) $(CFLAGS) $(FLAGS) $(opt) $(COMP) $(INCLUDE) -I$(BOOST_ROOT) -c $< -o $@
$(lib)/:
	mkdir -p $(lib)/
libSkin++.a: $(OBJC) $(lib)/
	ar r $@ $(OBJC)
	cp $@ $(lib)/$@
$(bin)/:
	mkdir -p $(bin)/
Skin++: libSkin++.a Skin++.o $(bin)/ sub-make
	$(F90) $(opt) $(FOMP) Skin++.o $(lib)/libSkin++.a $(libDon)/libDonjon.a $(libDra)/libDragon.a $(libTri)/libTrivac.a \
	$(libDra)/libDragon.a $(libUtl)/libUtilib.a $(libGan)/libGanlib.a $(clib) $(LFLAGS) -o Skin++
	cp $@ $(bin)/$@
clean:
	$(MAKE) -C ../../Donjon/src clean
	/bin/rm -f *.o *.a sub-make temp.* Skin++
