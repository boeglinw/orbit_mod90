#------Makefile for modules and the flux library ----------

# modules
MODOBJ = limiter_mod.o limiter_control_mod.o

INCLUDE = ../include/ -I../flux/

FLUX       = ../flux/
SPLINE     = ../spline/
RDEQDSK    = ../read_eqdsk/
MODULES    = ../include/
LIMITER    = ../limiter/

F90     = gfortran

#FFLAGS  = -ffixed-line-length-none -fd-lines-as-comments \

# FDFLAG  = -fd-lines-as-code

FFLAGS  = -ffixed-line-length-none $(FDFLAG) -fPIC \
          -w -fno-automatic -I$(INCLUDE)


THIS_LIB = liblimiter.a


#----------------------------------------------------------
all: $(MODOBJ) $(THIS_LIB)

%.o: %.f90
	$(F90) $(FFLAGS) -c $<


test_limiter: test_limiter.o $(MODOBJ)
	$(F90)  test_limiter.o  -L$(LIMITER) -llimiter -L$(FLUX) -lflux  -L$(RDEQDSK) -lread_eqdsk -L$(MODULES) -lmodules -L$(SPLINE) -lspline  -o test_limiter


$(THIS_LIB): $(MODOBJ)
	ar rv $(THIS_LIB) $(MODOBJ)
	ranlib $(THIS_LIB)


#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f *.o *.a *.mod test_limiter 




