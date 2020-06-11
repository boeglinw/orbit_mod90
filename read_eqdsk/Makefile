#------Makefile for rdeqdsk split into subroutines----------

# modules
# MODOBJ = ../include/orbit_parameters.o ../flux/flux_module.o 

RDEQDSK_OBJ = read_eqdsk.o
INCLUDE = ../include/ -I../flux/

F90     = gfortran

#FFLAGS  = -ffixed-line-length-none -fd-lines-as-comments \

# FDFLAG  = -fd-lines-as-code

FFLAGS  = -ffixed-line-length-none $(FDFLAG) \
          -w -fno-automatic -I$(INCLUDE)

SPLINE     = ../spline/
MODULES    = ../include/
FLUX       = ../flux/

THIS_LIB = libread_eqdsk.a


#----------------------------------------------------------
all: $(RDEQDSK_OBJ) $(THIS_LIB)

try: try.o $(RDEQDSK_OBJ)
	$(F90) -L$(MODULES) -lmodules -L$(FLUX) -lflux -L$(SPLINE) -lspline -L./ -lread_eqdsk try.o -o try

%.o: %.f90
	$(F90) -c $(FFLAGS) $< -o $@


$(THIS_LIB): $(RDEQDSK_OBJ)
	ar rv $(THIS_LIB) $(RDEQDSK_OBJ)
	ranlib $(THIS_LIB)

#example of creating object files the long way
#magfld.o: magfld.f
#	$(F90) -c magfld.f $(FFLAGS) -o $@
#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f $(RDEQDSK_OBJ)
	rm -f $(THIS_LIB)
	rm *.mod


