#FC = ifort
#FC = gfortran
FFLAGS = -g -C

default: hexa_test.exe matrix4d_test.exe

hexa_test.exe: hexa_test.o module_hexa.o module_io_matlab.o Makefile
	$(FC) $(FFLAGS) -o hexa_test.exe hexa_test.o module_hexa.o module_io_matlab.o

module_hexa.o: module_hexa.f90 Makefile
	$(FC) $(FFLAGS) -c module_hexa.f90

hexa_test.o: hexa_test.f90 module_hexa.o module_io_matlab.o Makefile
	$(FC) $(FFLAGS) -c hexa_test.f90

module_io_matlab.o: module_io_matlab.f90  Makefile
	$(FC) $(FFLAGS) -c module_io_matlab.f90

module_matrix4d.o: module_matrix4d.f90 Makefile
	$(FC) $(FFLAGS) -c module_matrix4d.f90

matrix4d_test.exe: matrix4d_test.f90 module_matrix4d.o module_io_matlab.o Makefile
	$(FC) $(FFLAGS) -o matrix4d_test.exe matrix4d_test.f90 module_matrix4d.o module_io_matlab.o


clean:
	rm -rf *.mod *.o *.exe
