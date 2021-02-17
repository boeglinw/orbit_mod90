#------Makefile for the boris_cylinder library ----------

# modules
MODOBJ = boris_cylinder.o

MODULES    = ../include/

INCLUDE = $(MODULES)

F90     = gfortran

FFLAGS  = -ffixed-line-length-none $(FDFLAG) -fPIC \
          -w -fno-automatic -I$(INCLUDE)

THIS_LIB = libboris_cylinder.a

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




