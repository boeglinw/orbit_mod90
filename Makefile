#-----Master Makefile for orbit_mod90-------

# module based lorentz orbit tracker

INCLUDE    = ./include/
SPLINE     = ./spline/
FLUX       = ./flux/
RDEQDSK    = ./read_eqdsk/
BORIS      = ./boris_pusher/
LIMITER      = ./limiter/
PYTHON      = ./python/
BSINTEGRATOR = ./bs_integrator/
TRACKER = ./tracker/


#----------------------------------------------------------
# targets need to be a different names that directory names
ALL:  include_mod spline_mod flux_mod  read_eqdsk_mod limiter_module boris  bs_module tracker_module

include_mod:
	cd $(INCLUDE) ; make all

spline_mod:
	cd $(SPLINE); make all

flux_mod:
	cd $(FLUX); make all

read_eqdsk_mod:
	cd $(RDEQDSK); make all

limiter_module:
	cd $(LIMITER); make all

boris:
	cd $(BORIS); make all

bs_module:
	cd $(BSINTEGRATOR); make all

tracker_module:
	cd $(TRACKER); make all



python_modules:
	cd $(FLUX); make -f flux.mak
	cd $(TRACKER); make -f Trpy.mak
	cd $(BORIS); make -f BorisCylpy.mak
	cp $(BORIS)/*.so $(PYTHON)/.
	cp $(FLUX)/*.so $(PYTHON)/.
	cp $(TRACKER)/*.so $(PYTHON)/.

#----------------------------------------------------------

clean: clean_include clean_spline clean_flux clean_read_eqdsk clean_boris clean_limiter clean_bs clean_tracker clean_python

clean_include:
	cd $(INCLUDE); make clean

clean_spline:
	cd $(SPLINE); make clean

clean_flux:
	cd $(FLUX); make clean; make -f flux.mak clean

clean_read_eqdsk:
	cd $(RDEQDSK); make clean

clean_boris:
	cd $(BORIS); make clean;  make -f BorisCylpy.mak clean

clean_limiter:
	cd $(LIMITER); make clean

clean_bs:
	cd $(BSINTEGRATOR); make clean

clean_tracker:
	cd $(TRACKER); make clean

clean_python:
	cd $(PYTHON); rm *.so *.npz
