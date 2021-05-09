FC = ifort
FC = gfortran
FFLAGS = -g -C --ffpe-trap=invalid,zero,overflow,underflow -traceback -fpe0 # intel debug
FFLAGS = -O3 
FFLAGS = -O3 -fdefault-real-8  # gnu
FFLAGS = -g -C -fpe0 -r8  # intel debug 
FFLAGS = -g -C  -fdefault-real-8

default: hexa_test.exe\
 ndt_mult_test.exe\
 sweeps_test.exe\
 f_assembly_test.exe\
 prolongation_test.exe\
 restriction_test.exe\
 ndt_assembly_test.exe\
 ndt_w_test.exe\
 coarsening_icl_test.exe\
 ndt_boundary_conditions_test.exe\
 vec_boundary_conditions_test.exe\
 femwind_test.exe



hexa_test.exe: hexa_test.o module_hexa.o module_utils.o Makefile
	$(FC) $(FFLAGS) -o hexa_test.exe hexa_test.o module_hexa.o module_utils.o

module_hexa.o: module_hexa.f90 module_utils.o Makefile
	$(FC) $(FFLAGS) -c module_hexa.f90

hexa_test.o: hexa_test.f90 module_hexa.o module_utils.o Makefile
	$(FC) $(FFLAGS) -c hexa_test.f90

module_utils.o: module_utils.f90  Makefile
	$(FC) $(FFLAGS) -c module_utils.f90

module_common.o: module_common.f90  module_utils.o Makefile
	$(FC) $(FFLAGS) -c module_common.f90

module_ndt_mult.o: module_ndt_mult.f90 Makefile
	$(FC) $(FFLAGS) -c module_ndt_mult.f90

ndt_mult_test.exe: ndt_mult_test.f90 module_ndt_mult.o module_utils.o Makefile
	$(FC) $(FFLAGS) -o ndt_mult_test.exe ndt_mult_test.f90 module_ndt_mult.o module_utils.o

sweeps_test.exe: sweeps_test.o module_sweeps.o module_utils.o Makefile
	$(FC) $(FFLAGS) -o sweeps_test.exe sweeps_test.f90 module_sweeps.o module_utils.o

module_sweeps.o: module_sweeps.f90 Makefile
	$(FC) $(FFLAGS) -c module_sweeps.f90

sweeps_test.o: sweeps_test.f90 module_sweeps.o module_utils.o Makefile
	$(FC) $(FFLAGS) -c sweeps_test.f90

ndt_assembly_test.exe: ndt_assembly_test.o module_ndt_assembly.o module_utils.o module_hexa.o  Makefile
	$(FC) $(FFLAGS) -o ndt_assembly_test.exe ndt_assembly_test.o  module_ndt_assembly.o module_utils.o module_hexa.o 

ndt_assembly_test.o: ndt_assembly_test.f90 module_hexa.o module_utils.o module_ndt_assembly.o Makefile
	$(FC) $(FFLAGS) -c ndt_assembly_test.f90 

module_ndt_assembly.o: module_ndt_assembly.f90 module_hexa.o Makefile
	$(FC) $(FFLAGS) -c module_ndt_assembly.f90  

f_assembly_test.exe: f_assembly_test.f90 module_f_assembly.o module_utils.o module_hexa.o  Makefile
	$(FC) $(FFLAGS) -o f_assembly_test.exe f_assembly_test.f90  module_f_assembly.o module_utils.o module_hexa.o

module_f_assembly.o: module_f_assembly.f90 module_hexa.o module_utils.o Makefile
	$(FC) $(FFLAGS) -c module_f_assembly.f90

module_coarsening.o: module_coarsening.f90 module_hexa.o Makefile
	$(FC) $(FFLAGS) -c module_coarsening.f90  

PROLONGATION = prolongation_test.f90  module_coarsening.o module_utils.o
prolongation_test.exe: $(PROLONGATION) 
	$(FC) $(FFLAGS) -o prolongation_test.exe $(PROLONGATION) 

RESTRICTION = restriction_test.f90  module_coarsening.o module_utils.o
restriction_test.exe: $(RESTRICTION) Makefile
	$(FC) $(FFLAGS) -o restriction_test.exe $(RESTRICTION) 

ICL = coarsening_icl_test.f90  module_coarsening.o module_utils.o
coarsening_icl_test.exe: $(ICL) Makefile
	$(FC) $(FFLAGS) -o $@ $(ICL) 

module_boundary_conditions.o: module_boundary_conditions.f90 module_utils.o Makefile
	$(FC) $(FFLAGS) -c module_boundary_conditions.f90  

NDT_BOUNDARY_CONDITIONS = ndt_boundary_conditions_test.f90  module_boundary_conditions.o module_utils.o
ndt_boundary_conditions_test.exe: $(NDT_BOUNDARY_CONDITIONS) Makefile
	$(FC) $(FFLAGS) -o $@ $(NDT_BOUNDARY_CONDITIONS) 

VEC_BOUNDARY_CONDITIONS = vec_boundary_conditions_test.f90  module_boundary_conditions.o module_utils.o
vec_boundary_conditions_test.exe: $(VEC_BOUNDARY_CONDITIONS) Makefile
	$(FC) $(FFLAGS) -o $@ $(VEC_BOUNDARY_CONDITIONS)

ndt_w_test.exe: ndt_w_test.o module_ndt_w_assembly.o module_utils.o module_hexa.o module_lin_alg.o Makefile
	$(FC) $(FFLAGS) -o ndt_w_test.exe ndt_w_test.o module_ndt_w_assembly.o module_utils.o module_hexa.o module_lin_alg.o

ndt_w_test.o: ndt_w_test.f90  module_hexa.o module_utils.o  module_ndt_w_assembly.o Makefile
	$(FC) $(FFLAGS) -c ndt_w_test.f90

module_ndt_w_assembly.o:  module_ndt_w_assembly.f90 module_hexa.o module_lin_alg.o Makefile
	$(FC) $(FFLAGS) -c module_ndt_w_assembly.f90

module_lin_alg.o:  module_lin_alg.f90 Makefile
	$(FC) $(FFLAGS) -c module_lin_alg.f90 

module_multigrid.o: module_multigrid.f90 module_ndt_assembly.o module_boundary_conditions.o module_common.o module_utils.o Makefile
	$(FC) $(FFLAGS) -c  module_multigrid.f90

module_femwind.o: module_femwind.f90 module_multigrid.o module_f_assembly.o module_ndt_w_assembly.o module_boundary_conditions.o module_common.o module_utils.o Makefile
	$(FC) $(FFLAGS) -c  module_femwind.f90

FEMWIND_TEST = femwind_test.f90\
  module_femwind.o\
  module_multigrid.o\
  module_ndt_assembly.o\
  module_ndt_w_assembly.o\
  module_ndt_mult.o\
  module_lin_alg.o\
  module_f_assembly.o\
  module_hexa.o\
  module_coarsening.o\
  module_sweeps.o\
  module_boundary_conditions.o\
  module_utils.o \
  module_common.o
femwind_test.exe: $(FEMWIND_TEST)
	$(FC) $(FFLAGS) -o $@ $(FEMWIND_TEST)

clean:
	rm -rf *.mod *.o *.exe *.dSYM
