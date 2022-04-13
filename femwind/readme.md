Matlab prototype
================

This is femwind, the 3rd generation mass consistent downscaling.
The code is at https://github.com/openwfm/wrf-fire-matlab/tree/quicwind/femwind

Description of the math of the finite element model is at
https://www.overleaf.com/read/ptfyhxgfnntn 
section Finite elements. There is only a brief description of the multigrid method.

The code in this directory is meant to be self-contained. Files
from the previous generations quicwind and quicwind/saddlepoint 
as well as utilities elsewhere in this repository are copied here as needed.

To run: in Matlab

close all         % get rid of open figure windows
p=femwind_test    % copy params to p, then modify p as desired
femwind_test(p)   % run with params settings from p
p=femwind_test(p)   % run with params settings from p

Jan Mandel, February 2021

Fortran implementation
======================

cd fortran
edit Makefile to set flags
make
cd ..
in Matlab:
femwind_run_fortran_rate_test  ! regression test

Jan Mandel, Angela Morrison, Evan Shapiro, May 2021

