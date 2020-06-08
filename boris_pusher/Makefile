#------Makefile for modules and the flux library ----------

# modules
MODOBJ = boris_mod.o BT.o

INCLUDE = ../include/ -I../flux/

FLUX       = ../flux/
SPLINE     = ../spline/
RDEQDSK    = ../read_eqdsk/
MODULES    = ../include/
BORIS    = ../boris_pusher/

F90     = gfortran

#FFLAGS  = -ffixed-line-length-none -fd-lines-as-comments \

# FDFLAG  = -fd-lines-as-code

FFLAGS  = -ffixed-line-length-none $(FDFLAG) \
          -w -fno-automatic -I$(INCLUDE)



THIS_LIB = libboris.a


#----------------------------------------------------------
all: $(MODOBJ) $(THIS_LIB)

%.o: %.f90
	$(F90) $(FFLAGS) -c $<


test_boris: test_boris.o $(MODOBJ)
	$(F90)  -L$(MODULES) -lmodules -L$(BORIS) -lboris -L$(FLUX) -lflux -L$(SPLINE) -lspline -L$(RDEQDSK) -lread_eqdsk  -L$(BORIS) -lboris test_boris.o -o test_boris

test_boris_c: test_boris_c.o $(MODOBJ)
	$(F90)  -L$(MODULES) -lmodules -L$(BORIS) -lboris -L$(FLUX) -lflux -L$(SPLINE) -lspline -L$(RDEQDSK) -lread_eqdsk  -L$(BORIS) -lboris test_boris_c.o -o test_boris_c

test_BT: test_BT.o $(MODOBJ)
	$(F90)  -L$(MODULES) -lmodules -L$(BORIS) -lboris -L$(FLUX) -lflux -L$(SPLINE) -lspline -L$(RDEQDSK) -lread_eqdsk  -L$(BORIS) -lboris test_BT.o -o test_BT


$(THIS_LIB): $(MODOBJ)
	ar rv $(THIS_LIB) $(MODOBJ)
	ranlib $(THIS_LIB)



#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f *.o *.a *.mod




