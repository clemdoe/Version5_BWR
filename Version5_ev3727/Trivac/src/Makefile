#---------------------------------------------------------------------------
#
#  Makefile for building the Trivac library and load module
#  Author : A. Hebert (2018-5-10)
#
#---------------------------------------------------------------------------
#
ARCH = $(shell uname -m)
ifneq (,$(filter $(ARCH),aarch64 arm64))
  nbit =
else
  ifneq (,$(filter $(ARCH),i386 i686))
    nbit = -m32
  else
    nbit = -m64
  endif
endif

DIRNAME = $(shell uname -sm | sed 's/[ ]/_/')
OS = $(shell uname -s | cut -d"_" -f1)
opt = -O -g
ifeq ($(openmp),1)
  COMP = -fopenmp
else
  COMP =
endif

ifeq ($(intel),1)
  fcompiler = ifort
  ccompiler = icc
else
  ifeq ($(nvidia),1)
    fcompiler = nvfortran
    ccompiler = nvc
  else
    ifeq ($(llvm),1)
      fcompiler = flang-new
      ccompiler = clang
    else
      fcompiler = gfortran
      ccompiler = gcc
    endif
  endif
endif

ifeq ($(OS),AIX)
  python_version_major := 2
else
  python_version_full := $(wordlist 2,4,$(subst ., ,$(shell python --version 2>&1)))
  python_version_major := $(word 1,${python_version_full})
  ifneq ($(python_version_major),2)
    python_version_major := 3
  endif
endif

ifeq ($(OS),Darwin)
  ifeq ($(openmp),1)
    ccompiler = gcc-14
  endif
  F90 = $(fcompiler)
  FFLAGS = $(nbit) -fPIC
  FFLAG77 = $(nbit) -fPIC
  LFLAGS = $(nbit)
else
ifeq ($(OS),Linux)
  F90 = $(fcompiler)
  FFLAGS = -Wall $(nbit) -fPIC
  FFLAG77 = -Wall $(nbit) -fPIC
  LFLAGS = $(nbit)
else
ifeq ($(OS),CYGWIN)
  F90 = $(fcompiler)
  FFLAGS = -Wall $(nbit) -fPIC
  FFLAG77 = -Wall $(nbit) -fPIC
  LFLAGS = $(nbit)
else
ifeq ($(OS),SunOS)
  fcompiler =
  MAKE = gmake
  F90 = f90
  FFLAGS = $(nbit) -s -ftrap=%none
  FFLAG77 = $(nbit) -s -ftrap=%none
  LFLAGS = $(nbit)
else
ifeq ($(OS),AIX)
  fcompiler =
  opt = -O4
  MAKE = gmake
  DIRNAME = AIX
  F90 = xlf90
  FFLAGS = -qstrict -qmaxmem=-1 -qsuffix=f=f90
  FFLAG77 = -qstrict -qmaxmem=-1 -qxlf77=leadzero -qfixed
  LFLAGS = -qstrict -bmaxdata:0x80000000 -qipa
else
  $(error $(OS) is not a valid OS)
endif
endif
endif
endif
endif
ifeq ($(fcompiler),gfortran)
  ifneq (,$(filter $(ARCH),i386 i686 x86_64))
    summary =
  else
    summary = -ffpe-summary=none
  endif
  ifeq ($(OS),Darwin)
    summary = -ffpe-summary=none
  endif
  FFLAGS += $(summary)
  FFLAG77 += -frecord-marker=4 $(summary)
endif

ifeq ($(intel),1)
  FFLAGS = -fPIC
  FFLAG77 = -fPIC
  lib = ../lib/$(DIRNAME)_intel
  libUtl = ../../Utilib/lib/$(DIRNAME)_intel
  libGan = ../../Ganlib/lib/$(DIRNAME)_intel
  bin = ../bin/$(DIRNAME)_intel
  INCLUDE = -I../../Ganlib/lib/$(DIRNAME)_intel/modules/
else
  ifeq ($(nvidia),1)
    lib = ../lib/$(DIRNAME)_nvidia
    libUtl = ../../Utilib/lib/$(DIRNAME)_nvidia
    libGan = ../../Ganlib/lib/$(DIRNAME)_nvidia
    bin = ../bin/$(DIRNAME)_nvidia
    INCLUDE = -I../../Ganlib/lib/$(DIRNAME)_nvidia/modules/
  else
    ifeq ($(llvm),1)
      lib = ../lib/$(DIRNAME)_llvm
      libUtl = ../../Utilib/lib/$(DIRNAME)_llvm
      libGan = ../../Ganlib/lib/$(DIRNAME)_llvm
      bin = ../bin/$(DIRNAME)_llvm
      INCLUDE = -I../../Ganlib/lib/$(DIRNAME)_llvm/modules/
      FFLAGS += -mmlir -fdynamic-heap-array
      LFLAGS += -lclang_rt.osx
    else
      lib = ../lib/$(DIRNAME)
      libUtl = ../../Utilib/lib/$(DIRNAME)
      libGan = ../../Ganlib/lib/$(DIRNAME)
      bin = ../bin/$(DIRNAME)
      INCLUDE = -I../../Ganlib/lib/$(DIRNAME)/modules/
    endif
  endif
endif

ifeq ($(hdf5),1)
  FLAGS += -DHDF5_LIB -I${HDF5_INC}
  FFLAGS += -I${HDF5_INC}
  LFLAGS += -L${HDF5_API} -lhdf5
endif

SRC77 = $(shell ls *.f)
ifeq ($(python_version_major),2)
  SRC90 = $(shell python ../../script/make_depend.py *.f90)
else
  SRC90 = $(shell python3 ../../script/make_depend_py3.py *.f90)
endif
OBJ90 = $(SRC90:.f90=.o)
OBJ77 = $(SRC77:.f=.o)
all : sub-make Trivac
ifeq ($(openmp),1)
	@echo 'Trivac: openmp is defined'
endif
ifeq ($(intel),1)
	@echo 'Trivac: intel is defined'
endif
ifeq ($(nvidia),1)
	@echo 'Trivac: nvidia is defined'
endif
ifeq ($(llvm),1)
	@echo 'Trivac: llvm is defined'
endif
ifeq ($(hdf5),1)
	@echo 'Trivac: hdf5 is defined'
endif
sub-make:
	$(MAKE) openmp=$(openmp) intel=$(intel) nvidia=$(nvidia) llvm=$(llvm) -C ../../Utilib/src
	$(MAKE) openmp=$(openmp) intel=$(intel) nvidia=$(nvidia) llvm=$(llvm) hdf5=$(hdf5) -C ../../Ganlib/src
%.o : %.f90
	$(F90) $(FFLAGS) $(opt) $(COMP) $(INCLUDE) -c $< -o $@
%.o : %.f
	$(F90) $(FFLAG77) $(opt) $(COMP) $(INCLUDE) -c $< -o $@
$(lib)/:
	mkdir -p $(lib)/
libTrivac.a: $(OBJ90) $(OBJ77) $(lib)/
	ar r $@ $(OBJ90) $(OBJ77)
	cp $@ $(lib)/$@
$(bin)/:
	mkdir -p $(bin)/
Trivac: libTrivac.a TRIVAC.o $(bin)/ sub-make
	$(F90) $(opt) $(COMP) TRIVAC.o $(lib)/libTrivac.a $(libUtl)/libUtilib.a \
	$(libGan)/libGanlib.a $(LFLAGS) -o Trivac
	cp $@ $(bin)/$@
clean:
	$(MAKE) -C ../../Utilib/src clean
	$(MAKE) -C ../../Ganlib/src clean
	/bin/rm -f *.o *.a sub-make temp.* Trivac
