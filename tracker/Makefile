#------Makefile for modules and the flux library ----------

# modules
MODOBJ = Tr.o 

FLUX       = ../flux/
SPLINE     = ../spline/
RDEQDSK    = ../read_eqdsk/
MODULES    = ../include/
LIMITER    = ../limiter/
BORIS    = ../boris_pusher/
BSINTEGRATOR = ../bs_integrator/
TRACKER = ../tracker/

INCLUDE = $(MODULES) -I $(FLUX) -I $(LIMITER) -I $(BORIS) -I $(BSINTEGRATOR) 

F90     = gfortran

#FFLAGS  = -ffixed-line-length-none -fd-lines-as-comments \

# FDFLAG  = -fd-lines-as-code

FFLAGS  = -ffixed-line-length-none -fPIC $(FDFLAG) \
          -w -fno-automatic -I$(INCLUDE)

THIS_LIB = libtracker.a

#----------------------------------------------------------
all: $(MODOBJ) $(THIS_LIB)

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

test_Tr: test_Tr.o $(MODOBJ)
	$(F90) test_Tr.o   -L$(TRACKER) -ltracker -L$(BORIS) -lboris  -L$(BSINTEGRATOR) -lbs_integrator -L$(LIMITER) -llimiter -L$(FLUX) -lflux  -L$(RDEQDSK) -lread_eqdsk -L$(MODULES) -lmodules -L$(SPLINE) -lspline -o test_Tr


$(THIS_LIB): $(MODOBJ)
	ar rv $(THIS_LIB) $(MODOBJ)
	ranlib $(THIS_LIB)



#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f *.o *.a *.mod test_Tr




