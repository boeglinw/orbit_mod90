#------Makefile for wraping the boris tracker library ----------

# modules
SRC = Trpy.f90
OBJ = Trpy.o

INCLUDE    = ../include/
BORIS      = ../boris_pusher/
FLUX       = ../flux/
SPLINE     = ../spline/
RDEQDSK    = ../read_eqdsk/
MODULES    = ../include/
LIMITER    = ../limiter/
BSINTEGRATOR = ../bs_integrator/
TRACKER = ../tracker/

F2PY     = f2py

F2PY_F1   = --include-paths $(INCLUDE):$(FLUX):$(LIMITER):$(BORIS):$(BSINTEGRATOR) --overwrite-signature -m 
F2PY_F2   = -c --fcompiler=gfortran --f90flags="-ffixed-line-length-none -w -fno-automatic -I../include/ -I../flux/ -I../limiter/ -I../boris_pusher/ -I../bs_integrator/ "

PROGRAM = Trpy

LIBS = -L$(TRACKER) -ltracker -L$(BORIS) -lboris  -L$(BSINTEGRATOR) -lbs_integrator -L$(LIMITER) -llimiter -L$(FLUX) -lflux  -L$(RDEQDSK) -lread_eqdsk -L$(MODULES) -lmodules -L$(SPLINE) -lspline

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
	rm -f *.so *.pyf fluxpy.o
	rm -rf *.dSYM




