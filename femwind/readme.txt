This is femwind, the 3rd generation mass consistent downscaling.
The code is at https://github.com/openwfm/wrf-fire-matlab/tree/quicwind/femwind

Description of the math of the finite element model is at
https://www.overleaf.com/read/ptfyhxgfnntn 
section Finite elements. There is no description of the multigrid method yet.
The method implemented now is two-level multigrid only.

The code in this directory is meant to be self-contained. Files
from the previous generations quicwind and quicwind/saddlepoint 
as well as utilities elsewhere in this repository are copied here as needed.

To run: in Matlab

close all
clear all
femwind_test
