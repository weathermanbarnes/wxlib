
LIB=../lib/libspherepack.a

UNAMES := $(shell uname -s)

ifeq ($(UNAMES),Linux)

  F90 := gfortran -fdefault-real-8 -fPIC -fallow-argument-mismatch -I../lib
#  F90 := gfortran -fdefault-real-8 -ffpe-trap=invalid,zero,overflow -fbacktrace -g -fPIC -I../lib
  CPP := gfortran -cpp
  MAKE := make
  AR := /usr/bin/ar

endif

ifeq ($(UNAMES),AIX)

  F90 := xlf -qmoddir=../lib -I../lib
  CPP := xlf -d -qnoobject
  MAKE := gmake
  AR := /usr/bin/ar

endif

ifeq ($(UNAMES),SunOS)

    AR := /usr/ccs/bin/ar
    F90 := /opt/SUNWspro/bin/f90 -moddir=../lib -I../lib
    CPP := /opt/SUNWspro/bin/f90 -F
    MAKE := /fs/local/bin/make

endif

ifeq ($(UNAMES),IRIX64)

    AR := /usr/bin/ar
    F90 := f90 -I../lib
    CPP := f90 -E
    MAKE := /usr/local/bin/gmake

endif

ifeq ($(UNAMES),Darwin)

    AR := /usr/bin/ar
    F90 := gfortran 
#   F90 := gfortran -fdefault-real-8 
    CPP := gfortran -cpp
    MAKE := /usr/bin/gmake

endif
