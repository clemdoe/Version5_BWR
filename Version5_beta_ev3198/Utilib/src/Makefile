#---------------------------------------------------------------------------
#
#  Makefile for building the Utilib library
#  Author : A. Hebert (2018-5-10)
#
#---------------------------------------------------------------------------
#
ARCH = $(shell uname -m)
ifeq ($(ARCH),aarch64)
  nbit =
  summary = -ffpe-summary=none
else
  ifeq (,$(filter $(ARCH),i386 i686))
    nbit = -m64
    summary = -ffpe-summary=none
  else
    nbit = -m32
    summary =
  endif
endif

DIRNAME = $(shell uname -sm | sed 's/[ ]/_/')
OS = $(shell uname -s | cut -d"_" -f1)
ccompiler = gcc
opt = -O -g
ifeq ($(openmp),1)
  COMP = -fopenmp
else
  COMP =
endif

ifeq ($(intel),1)
  fcompiler = ifort
else
  fcompiler = gfortran
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
  ccompiler = gcc-12
endif
  F90 = $(fcompiler)
  C = $(ccompiler)
  CFLAGS = -Wall $(nbit) -fPIC
  FFLAGS = -Wall $(nbit) -fPIC $(summary)
  FFLAG77 = -Wall $(nbit) -fPIC $(summary) -frecord-marker=4
else
ifeq ($(OS),Linux)
  F90 = $(fcompiler)
  C = $(ccompiler)
  CFLAGS = -Wall $(nbit) -fPIC
  FFLAGS = -Wall $(nbit) -fPIC
  FFLAG77 = -Wall $(nbit) -fPIC -frecord-marker=4
else
ifeq ($(OS),CYGWIN)
  F90 = $(fcompiler)
  C = $(ccompiler)
  CFLAGS = -Wall $(nbit) -fPIC
  FFLAGS = -Wall $(nbit) -fPIC
  FFLAG77 = -Wall $(nbit) -fPIC -frecord-marker=4
else
ifeq ($(OS),SunOS)
  F90 = f90
  C = cc
  CFLAGS = $(nbit)
  FFLAGS = $(nbit) -s -ftrap=%none
  FFLAG77 = $(nbit) -s -ftrap=%none
else
ifeq ($(OS),AIX)
  opt = -O4
  DIRNAME = AIX
  F90 = xlf90
  C = cc
  CFLAGS = -qstrict
  FFLAGS = -qstrict -qmaxmem=-1 -qsuffix=f=f90
  FFLAG77 = -qstrict -qmaxmem=-1 -qxlf77=leadzero -qfixed
else
  $(error $(OS) is not a valid OS)
endif
endif
endif
endif
endif

ifeq ($(intel),1)
  FFLAGS = -fPIC
  FFLAG77 = -fPIC
  lib = ../lib/$(DIRNAME)_intel
  lib_module = ../lib/$(DIRNAME)_intel/modules
else
  lib = ../lib/$(DIRNAME)
  lib_module = ../lib/$(DIRNAME)/modules
endif

SRCC = $(shell ls *.c)
SRC77 = $(shell ls *.f)
ifeq ($(python_version_major),2)
  SRC90 = $(shell python ../../script/make_depend.py *.f90)
else
  SRC90 = $(shell python3 ../../script/make_depend_py3.py *.f90)
endif
OBJC = $(SRCC:.c=.o)
OBJ90 = $(SRC90:.f90=.o)
OBJ77 = $(SRC77:.f=.o)
all : libUtilib.a
ifeq ($(openmp),1)
	@echo 'Utilib: openmp is defined'
endif
ifeq ($(intel),1)
	@echo 'Utilib: intel is defined'
endif
	@echo "Utilib: python version=" ${python_version_major}
%.o : %.c
	$(C) $(CFLAGS) $(opt) $(COMP) -c $< -o $@
%.o : %.f90
	$(F90) $(FFLAGS) $(opt) $(COMP) $(INCLUDE) -c $< -o $@
%.o : %.f
	$(F90) $(FFLAG77) $(opt) $(COMP) -c $< -o $@
all: $(lib_module)
$(lib_module)/:
	mkdir -p $(lib_module)/
$(lib)/: $(lib_module)/
	mkdir -p $(lib)/
libUtilib.a: $(OBJC) $(OBJ90) $(OBJ77) $(lib)/
	ar r $@ $(OBJC) $(OBJ90) $(OBJ77)
	cp $@ $(lib)/$@
clean:
	/bin/rm -f *.o *.mod *.a
