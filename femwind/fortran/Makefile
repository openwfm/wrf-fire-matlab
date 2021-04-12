#FC = ifort
#FC = gfortran
#FFLAGS = -g -C --ffpe-trap=invalid,zero,overflow,underflow
FFLAGS = -g -C -fpe0 -traceback

default: hexa_test.exe ndt_mult_test.exe  sweeps_test.exe ndt_assembly_test.exe prolongation_test.exe restriction_test.exe

hexa_test.exe: hexa_test.o module_hexa.o module_io_matlab.o Makefile
	$(FC) $(FFLAGS) -o hexa_test.exe hexa_test.o module_hexa.o module_io_matlab.o

module_hexa.o: module_hexa.f90 Makefile
	$(FC) $(FFLAGS) -c module_hexa.f90

hexa_test.o: hexa_test.f90 module_hexa.o module_io_matlab.o Makefile
	$(FC) $(FFLAGS) -c hexa_test.f90

module_io_matlab.o: module_io_matlab.f90  Makefile
	$(FC) $(FFLAGS) -c module_io_matlab.f90

module_ndt_mult.o: module_ndt_mult.f90 Makefile
	$(FC) $(FFLAGS) -c module_ndt_mult.f90

ndt_mult_test.exe: ndt_mult_test.f90 module_ndt_mult.o module_io_matlab.o Makefile
	$(FC) $(FFLAGS) -o ndt_mult_test.exe ndt_mult_test.f90 module_ndt_mult.o module_io_matlab.o

sweeps_test.exe: sweeps_test.o module_sweeps.o module_io_matlab.o Makefile
	$(FC) $(FFLAGS) -o sweeps_test.exe sweeps_test.f90 module_sweeps.o module_io_matlab.o

module_sweeps.o: module_sweeps.f90 Makefile
	$(FC) $(FFLAGS) -c module_sweeps.f90

sweeps_test.o: sweeps_test.f90 module_sweeps.o module_io_matlab.o Makefile
	$(FC) $(FFLAGS) -c sweeps_test.f90

ndt_assembly_test.exe: ndt_assembly_test.o module_ndt_assembly.o module_io_matlab.o module_hexa.o  Makefile
	$(FC) $(FFLAGS) -o ndt_assembly_test.exe ndt_assembly_test.o  module_ndt_assembly.o module_io_matlab.o module_hexa.o 

ndt_assembly_test.o: ndt_assembly_test.f90 module_hexa.o module_io_matlab.o module_ndt_assembly.o Makefile
	$(FC) $(FFLAGS) -c ndt_assembly_test.f90 

module_ndt_assembly.o: module_ndt_assembly.f90 module_hexa.o Makefile
	$(FC) $(FFLAGS) -c module_ndt_assembly.f90  

module_coarsening.o: module_coarsening.f90 module_hexa.o Makefile
	$(FC) $(FFLAGS) -c module_coarsening.f90  

PROLONGATION = prolongation_test.f90  module_coarsening.o module_io_matlab.o
prolongation_test.exe: $(PROLONGATION) 
	$(FC) $(FFLAGS) -o prolongation_test.exe $(PROLONGATION) 

RESTRICTION = restriction_test.f90  module_coarsening.o module_io_matlab.o
restriction_test.exe: $(RESTRICTION) 
	$(FC) $(FFLAGS) -o restriction_test.exe $(RESTRICTION) 

clean:
	rm -rf *.mod *.o *.exe
