#------Makefile for modules and the flux library ----------

# modules
MODOBJ = flux_par_mod.o em_fields_mod.o bfield.o checklim.o get_flux.o is_inside.o get_psi.o

INCLUDE = ../include/ -I../flux/

FLUX       = ../flux/
SPLINE     = ../spline/
RDEQDSK    = ../read_eqdsk/
MODULES    = ../include/

F90     = gfortran

#FFLAGS  = -ffixed-line-length-none -fd-lines-as-comments \

# FDFLAG  = -fd-lines-as-code



FFLAGS  = -ffixed-line-length-none -fPIC $(FDFLAG) \
          -w -fno-automatic -I$(INCLUDE)



THIS_LIB = libflux.a


#----------------------------------------------------------
all: $(MODOBJ) $(THIS_LIB)

%.o: %.f90
	$(F90) $(FFLAGS) -c $<


test_bfield: test_bfield.o $(MODOBJ)
	$(F90)  test_bfield.o -L$(FLUX) -lflux -L$(RDEQDSK) -lread_eqdsk  -L$(MODULES) -lmodules -L$(SPLINE) -lspline  -o test_bfield


test_checklim: test_checklim.o $(MODOBJ)
	$(F90)  test_checklim.o  -L$(FLUX) -lflux  -L$(RDEQDSK) -lread_eqdsk -L$(MODULES) -lmodules -L$(SPLINE) -lspline -o test_checklim


$(THIS_LIB): $(MODOBJ)
	ar rv $(THIS_LIB) $(MODOBJ)
	ranlib $(THIS_LIB)



#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f *.o *.a *.mod test_bfield test_zlim




