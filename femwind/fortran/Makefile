#FC = gfortran
#CC = gcc
#FC = ifort
#CC = icc
FFLAGS = -g -C --ffpe-trap=invalid,zero,overflow,underflow -traceback -fpe0 # intel debug
FFLAGS = -O3 -fdefault-real-8  # gnu
FFLAGS = -g -C -fpe0 -r8  # intel debug 
FFLAGS = -g -C  -fdefault-real-8
FFLAGS = -O3 -pg
FFLAGS = -g -C -I /opt/local/include
LFLAGS = -L/opt/local/lib -lnetcdf -l netcdff
FFLAGS = -pg -O3 -I $(NETCDF)/include
FFLAGSOPT = -O3 -I $(NETCDF)/include
FFLAGS = -g -C -I $(NETCDF)/include 
FFLAGS = -O3 -I $(NETCDF)/include
LFLAGS = -L$(NETCDF)/lib -lnetcdf -lnetcdff

default: all

wrf: femwind_wrfout.exe

tests: clean_txt\
 hexa_test.exe\
 ndt_mult_test.exe\
 sweeps_test.exe\
 f_assembly_test.exe\
 prolongation_test.exe\
 restriction_test.exe\
 ndt_assembly_test.exe\
 w_assembly_test.exe\
 coarsening_icl_test.exe\
 ndt_boundary_conditions_test.exe\
 vec_boundary_conditions_test.exe\
 femwind_test.exe\
 netcdf_test.exe\
 msleep_test.exe

all: wrf tests

hexa_test.exe: hexa_test.o module_hexa.o module_utils.o Makefile
	$(FC) $(FFLAGS) -o hexa_test.exe hexa_test.o module_hexa.o module_utils.o

module_hexa.o: module_hexa.f90 module_utils.o Makefile
	$(FC) $(FFLAGSOPT) -c module_hexa.f90

hexa_test.o: hexa_test.f90 module_hexa.o module_utils.o Makefile
	$(FC) $(FFLAGS) -c hexa_test.f90

module_utils.o: module_utils.f90  Makefile
	$(FC) $(FFLAGS) -c module_utils.f90

module_netcdf.o: module_utils.o module_netcdf.f90 Makefile
	$(FC) $(FFLAGS) -c module_netcdf.f90

module_wrfout.o: module_wrfout.f90 module_netcdf.o module_msleep.o Makefile
	$(FC) $(FFLAGS) -c module_wrfout.f90

module_common.o: module_common.f90  module_utils.o Makefile
	$(FC) $(FFLAGS) -c module_common.f90

module_ndt_mult.o: module_ndt_mult.f90  module_utils.o Makefile
	$(FC) $(FFLAGS) -c module_ndt_mult.f90

ndt_mult_test.exe: ndt_mult_test.f90 module_ndt_mult.o module_utils.o Makefile
	$(FC) $(FFLAGS) -o ndt_mult_test.exe ndt_mult_test.f90 module_ndt_mult.o module_utils.o

sweeps_test.exe: sweeps_test.o module_sweeps.o module_utils.o Makefile
	$(FC) $(FFLAGS) -o sweeps_test.exe sweeps_test.f90 module_sweeps.o module_utils.o

module_sweeps.o: module_sweeps.f90 module_utils.o Makefile
	$(FC) $(FFLAGS) -c module_sweeps.f90

sweeps_test.o: sweeps_test.f90 module_sweeps.o module_utils.o Makefile
	$(FC) $(FFLAGS) -c sweeps_test.f90

ndt_assembly_test.exe: ndt_assembly_test.o module_ndt_assembly.o module_utils.o module_hexa.o  Makefile
	$(FC) $(FFLAG) -o ndt_assembly_test.exe ndt_assembly_test.o  module_ndt_assembly.o module_utils.o module_hexa.o 

ndt_assembly_test.o: ndt_assembly_test.f90 module_hexa.o module_utils.o module_ndt_assembly.o Makefile
	$(FC) $(FFLAGS) -c ndt_assembly_test.f90 

module_ndt_assembly.o: module_ndt_assembly.f90 module_hexa.o Makefile
	$(FC) $(FFLAGSOPT) -c module_ndt_assembly.f90  

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

w_assembly_test.exe: w_assembly_test.o module_w_assembly.o module_utils.o module_hexa.o module_lin_alg.o Makefile
	$(FC) $(FFLAGS) -o w_assembly_test.exe w_assembly_test.o module_w_assembly.o module_utils.o module_hexa.o module_lin_alg.o

w_assembly_test.o: w_assembly_test.f90  module_hexa.o module_utils.o  module_w_assembly.o Makefile
	$(FC) $(FFLAGS) -c w_assembly_test.f90

module_w_assembly.o:  module_w_assembly.f90 module_hexa.o module_lin_alg.o Makefile
	$(FC) $(FFLAGS) -c module_w_assembly.f90

module_lin_alg.o:  module_lin_alg.f90 Makefile
	$(FC) $(FFLAGS) -c module_lin_alg.f90 

module_multigrid.o: module_sweeps.o module_multigrid.f90 module_coarsening.o module_ndt_mult.o module_ndt_assembly.o module_boundary_conditions.o module_common.o module_utils.o Makefile
	$(FC) $(FFLAGS) -c  module_multigrid.f90

module_femwind.o: module_femwind.f90 module_multigrid.o module_f_assembly.o module_w_assembly.o module_boundary_conditions.o module_common.o module_utils.o Makefile
	$(FC) $(FFLAGS) -c  module_femwind.f90

FEMWIND_MODULES = \
  module_femwind.o\
  module_multigrid.o\
  module_ndt_assembly.o\
  module_w_assembly.o\
  module_ndt_mult.o\
  module_lin_alg.o\
  module_f_assembly.o\
  module_hexa.o\
  module_coarsening.o\
  module_sweeps.o\
  module_boundary_conditions.o\
  module_utils.o \
  module_common.o

femwind_test.exe: femwind_test.f90 $(FEMWIND_MODULES) Makefile
	$(FC) $(FFLAGS) -o $@ femwind_test.f90 $(FEMWIND_MODULES) 

WRFOUT_MODULES = \
  $(MSLEEP) \
  module_netcdf.o \
  module_wrfout.o

FEMWIND_WRFOUT = femwind_wrfout.f90 $(WRFOUT_MODULES) $(FEMWIND_MODULES)

femwind_wrfout.exe: $(FEMWIND_WRFOUT) Makefile
	$(FC) $(FFLAGS) -o $@ $(FEMWIND_WRFOUT) $(LFLAGS)

NETCDF_TEST = netcdf_test.f90 $(WRFOUT_MODULES) module_utils.o 

netcdf_test.exe: $(NETCDF_TEST) Makefile
	$(FC) $(FFLAGS) -o $@ $(NETCDF_TEST) $(LFLAGS)

msleep.o: msleep.c
	$(CC) -c msleep.c 

module_msleep.o: module_msleep.f90
	$(FC) $(FFLAGS) -c module_msleep.f90

MSLEEP = \
  module_msleep.o \
  msleep.o

msleep_test.exe: msleep_test.f90 $(MSLEEP)
	$(FC) -o msleep_test.exe msleep_test.f90 $(MSLEEP)

# delete test data files 
clean_txt:
	rm -f ../*.txt

clean: clean_txt
	rm -f *.mod *.o *.exe *.dSYM 
