#------Makefile for modules ----------

# modules


MODOBJ = orbit_parameters_mod.o orbit_detectors_mod.o  constants_and_masses_mod.o helper_functions_mod.o control_mod.o

INCLUDE = ../include/ -I../flux/

F90     = gfortran

#FFLAGS  = -ffixed-line-length-none -fd-lines-as-comments \

# FDFLAG  = -fd-lines-as-code

FFLAGS  = -ffixed-line-length-none -fPIC $(FDFLAG) \
          -w -fno-automatic -I$(INCLUDE)

THIS_LIB = libmodules.a

#----------------------------------------------------------
all: $(MODOBJ) $(THIS_LIB)

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

$(THIS_LIB): $(MODOBJ)
	ar rv $(THIS_LIB) $(MODOBJ)
	ranlib $(THIS_LIB)

#example of creating object files the long way
#magfld.o: magfld.f
#	$(F90) -c magfld.f $(FFLAGS) -o $@
#----------------------------------------------------------
.PHONY : clean

orbit_detectors_mod.o: orbit_detectors_mod.f90 
	$(F90) $(FFLAGS) -c $<
    
clean:
	rm -f $(MODOBJ) $(THIS_LIB) *.mod


