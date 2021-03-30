#FC = ifort
#FC = gfortran
#FFLAGS = -g -C --ffpe-trap=invalid,zero,overflow,underflow
FFLAGS = -g -C -fpe0 -traceback

default: hexa_test.exe ndt_mult_test.exe

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

ndt_assembly.exe: ndt_assembly.o module_hexa.o  Makefile
    $(FC) $(FFLAGS) -o ndt_assembly.exe module_hexa.o 

ndt_assembly.o: ndt_assembly.f90 module_hexa.o Makefile
    $(FC) $(FFLAGS) -c ndt_assembly.f90


clean:
	rm -rf *.mod *.o *.exe
