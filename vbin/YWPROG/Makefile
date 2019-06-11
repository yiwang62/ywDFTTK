#
.SUFFIXES: .inc .f .f90 .F .c .cc
FC=pgf90
CC = g++
# fortran linker
FCL=$(FC)

#IMSL = -L /usr/global/imsl-5.00/CTT4.0/lib/lib.linux -limsl
IMSL = -L/usr/global/vni/CTT6.0/lib/lib.pgi64 -limsl

#FFLAGS = -O2  -tp p6 $(IMSL)
FFLAGS = -O2 -r8 $(F90FLAGS) $(LINK_F90) $(IMSL)

#set ALL `ls *.f | sed 's/\.f//g' | echo`

ALL = atafre atatCLD atatget atatsum atatCCC atatdos atatPrs pcoor


EXE = ${ALL}

all: $(EXE)
