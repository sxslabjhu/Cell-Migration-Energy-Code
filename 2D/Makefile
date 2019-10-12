FComp = ifort
#FComp = gfortran

#DEBUG = yes
DEBUG = no

CC  = gcc
C++ = g++
#=========================================================================


FILES=parameters.f90  myfft.f90 linsys.f90 geometry.f90 chemical_mod.f90 \
	solver_mod.f90 network_mod.f90 IBforce.f90 IBmod.f90 energy.f90 \
	mainimp.f90 Makefile input.par

OBJN=parameters.o  myfft.o geometry.o solver_mod.o IBforce.o IBmod.o \
	linsys.o chemical_mod.o network_mod.o energy.o mainimp.o

ifeq ($(FComp),ifort)
  FC     = $(FComp)
  ##FFLAG  =-O2 -fp-model strict -assume protect_parens -check bound
  ifeq ($(DEBUG),yes)
    FFLAG  =-g -r8 -assume protect_parens 
    ##FFLAG  =-g -r8 -fp-model strict -assume protect_parens -check bound
  else
    FFLAG  =-O2 -r8 -fp-model strict -assume protect_parens -check bound
  endif
  FFTWINC=-I$(HOME)/local/fftw/ifc/include
  FFTWDIR=-L$(HOME)/local/fftw/ifc/lib
  LIBDIR =-L$(HOME)/lib
  ##LIBS   =-lfftpack5.1_ifc -llapack_ifc -lblas_ifc
  LIBS   = -llapack_ifc -lblas_ifc
  CFLAG  =-O2
  CPPFLAG=-O2
endif
ifeq ($(FComp),gfortran)
  FC     = $(FComp)
  ifeq ($(DEBUG),yes)
    ##FFLAG  =-g -frecursive -fcheck=bounds,mem,pointer  -fdefault-real-8 -fdefault-double-8 -freal-4-real-8
    FFLAG  =-g -fcheck=bounds,mem,pointer  -fdefault-real-8 -fdefault-double-8 -freal-4-real-8
#FFLAG=-O2 -g -frecursive
  else
    FFLAG  =-O2 -frecursive -fcheck=bounds,mem,pointer  -fdefault-real-8 -fdefault-double-8 -freal-4-real-8
  endif
  FFTWINC=-I${HOME}/local/fftw/gfc/include/
  FFTWDIR=-L${HOME}/local/fftw/gfc/lib/
  LIBDIR =-L${HOME}/lib
  ##LIBS   =-lfftpack5.1 -llapack_gfc -lblas_gfc
  LIBS   = -llapack_gfc -lblas_gfc
  CFLAG  =-O2
  CPPFLAG=-O2
endif


.SUFFIXES: .f90 .f .c

default: pack.zip
#default: imp
all: clean imp
##all: clean imp pc pack.zip  sclean


pack.zip: $(FILES)
	-@cp pack.zip pack.zip.old
	@zip -9 -u -y pack.zip $(FILES)

imp: $(OBJN)
	$(FC) $(FFLAG) -o $@  $(OBJN) $(FFTWINC) $(FFTWDIR) -lfftw3  $(LIBDIR) $(LIBS)
	
	

clean:
	@rm -f $(OBJS) imp pc *.dat dump.m *.mod *.o
	@rm -fr imp.dSYM pc.dSYM ns.dSYM
sclean:
	@rm -f $(OBJS) imp pc *.dat dump.m *.mod *.o pack.zip
	@rm -fr imp.dSYM pc.dSYM ns.dSYM

%.mod : %.o
	@if [! -f $@ ]; then \
	  rm $< \
	  $(FC) $(FFLAG) -c -o $@ $<
	fi


%.o : %.f90
	@echo "Compiling $<"
	$(FC) $(FFLAG) -c -o $@ $<


.c.o:
	$(CC) -c $(CFLAG) $*.c
.C.o:
	$(C++) -c $(CPPFLAG) $*.C
