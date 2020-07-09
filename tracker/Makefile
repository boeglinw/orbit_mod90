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

FFLAGS  = -ffixed-line-length-none $(FDFLAG) \
          -w -fno-automatic -I$(INCLUDE)

THIS_LIB = libtracker.a

#----------------------------------------------------------
all: $(MODOBJ) $(THIS_LIB)

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

test_Tr: test_Tr.o $(MODOBJ)
	$(F90)  -L$(MODULES) -lmodules -L$(BORIS) -lboris -L$(FLUX) -lflux -L$(SPLINE) -lspline -L$(RDEQDSK) -lread_eqdsk  -L$(TRACKER) -ltracker -L$(BORIS) -lboris -L$(LIMITER) -llimiter -L$(BSINTEGRATOR) -lbs_integrator -L$(TRACKER) -ltracker test_Tr.o -o test_Tr


$(THIS_LIB): $(MODOBJ)
	ar rv $(THIS_LIB) $(MODOBJ)
	ranlib $(THIS_LIB)



#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f *.o *.a *.mod test_Tr




