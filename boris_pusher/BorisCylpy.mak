#------Makefile for the boris_cylinder Python Interface BorisCylpy ----------

# modules
SRC = BorisCylpy.f90
OBJ = BorisCylpy.o

INCLUDE    = ../include/

F2PY     = f2py

F2PY_F1   = --include-paths $(INCLUDE) --overwrite-signature -m 
F2PY_F2   = -c --fcompiler=gfortran --f90flags="-ffixed-line-length-none -w -fno-automatic  -I../include/"

PROGRAM = BorisCylpy

LIBS = -L$(INCLUDE) 

#----------------------------------------------------------
all: $(OBJ)

%.o: %.f90
	$(F2PY) $(F2PY_F1) $(PROGRAM) -h sgn_$(PROGRAM).pyf $<
	# need to modify one line in the signature file to make sure that the bfield_array works correctly, of the code changes this needs to be adjusted
	# cat sgn_fluxpy.pyf | sed s/"real(kind=8), allocatable,dimension(:,:) :: b_array"/"real(kind=8), allocatable,dimension(size(r),4) :: b_array"/ > xx
	# mv xx sgn_$(PROGRAM).pyf
	$(F2PY) $(F2PY_F2) sgn_$(PROGRAM).pyf $(SRC) $(LIBS)

#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f *.so *.pyf BorisCylpy.o
	rm -rf *.dSYM




