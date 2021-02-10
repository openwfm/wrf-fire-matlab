#FC = ifort
#FC = gfortran
FFLAGS = -g -C

default: hexa_test.exe
	./hexa_test.exe

hexa_test.exe: hexa_test.o module_hexa.o Makefile
	$(FC) $(FFLAGS) -o hexa_test.exe hexa_test.o module_hexa.o

module_hexa.o: module_hexa.f90 Makefile
	$(FC) $(FFLAGS) -c module_hexa.f90

hexa_test.o: hexa_test.f90 module_hexa.o Makefile
	$(FC) $(FFLAGS) -c hexa_test.f90

clean:
	rm -rf *.mod *.o *.exe
