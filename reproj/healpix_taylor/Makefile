HPXDIR = $(HEALPIX)
HPXINC = $(HPXDIR)/include
HPXLIB = $(HPXDIR)/lib

FC90    = ifort
CFLAGS  = -O3 -I. -cm -w -vec_report0 -sox -openmp -openmp_report0 -CB -traceback -fPIC
FFLAGS  = -lhealpix -lcfitsio -CB -traceback -O3 -fPIC

LIBLINK = -L$(HPXLIB) -L$(FITSIO)/lib -L$(MKLLIB) -mkl

#PRJ     = test
PRJ     = driver
OBJ     = healpix_taylor.o hpxtay_utils.o

default: $(PRJ)

$(PRJ): $(PRJ).o $(OBJ)
	$(FC90) -o hpx_taylor $(PRJ).o $(OBJ) $(LIBLINK) $(FFLAGS)

$(PRJ).o: $(PRJ).f90 $(OBJ)
	$(FC90) $(CFLAGS) -I$(HPXINC) -c $(PRJ).f90

%.o: %.F90
	$(FC90) $(CFLAGS) -I$(HPXINC) -c -o $*.o $*.F90

clean:
	rm -f *.mod *.o hpx_taylor test
