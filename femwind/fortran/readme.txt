1. log in to computer with intel fortran and modules, like cheyenne or kingspeak
2. module load intel
3. make

To test with gfortran only, change in the Makefile FC to gfortran

Future: FC should really be set in the environment not in the Makefile
which then requires no change to change compilers, the flags are close enough
