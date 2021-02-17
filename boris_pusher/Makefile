#------Makefile for modules and the flux library ----------

# modules
MODOBJ = boris_mod.o

FLUX       = ../flux/
SPLINE     = ../spline/
RDEQDSK    = ../read_eqdsk/
MODULES    = ../include/
LIMITER    = ../limiter/
BORIS    = ../boris_pusher/

INCLUDE = $(MODULES) -I $(FLUX) -I $(LIMITER)

F90     = gfortran

#FFLAGS  = -ffixed-line-length-none -fd-lines-as-comments \

# FDFLAG  = -fd-lines-as-code

FFLAGS  = -ffixed-line-length-none -fPIC $(FDFLAG) \
          -w -fno-automatic -I$(INCLUDE)

THIS_LIB = libboris.a

#----------------------------------------------------------
all: $(MODOBJ) $(THIS_LIB)

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

$(THIS_LIB): $(MODOBJ)
	ar rv $(THIS_LIB) $(MODOBJ)
	ranlib $(THIS_LIB)



#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f *.o *.a *.mod




