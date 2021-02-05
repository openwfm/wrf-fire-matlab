FC = ifort
FFLAGS = -g -C

hexa_test.exe: hexa_test.o module_hexa.o
	$(FC) $(FFLAGS) -o hexa_test.exe hexa_test.o module_hexa.o

module_hexa.o: module_hexa.f90
	$(FC) $(FFLAGS) -c module_hexa.f90

hexa_test.o: hexa_test.f90 module_hexa.o
	$(FC) $(FFLAGS) -c hexa_test.f90

clean:
	rm -rf *.mod *.o *.exe
