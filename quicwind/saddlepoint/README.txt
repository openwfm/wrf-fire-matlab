This is the 2nd generation quicwind code. It is currently hosted at
https://github.com/openwfm/wrf-fire-matlab/tree/quicwind/quicwind/saddlepoint
To see what it can do, try the current driver saddle_test_n_3d.m
The code in this directory is meant to be self-contained.
We are copying files from the parent directory as needed.
The description is still at https://www.overleaf.com/read/ptfyhxgfnntn

From Angel Dec 2020:

now I remember what tests we did:
saddle_test_3d.m: Old version of saddlepoint solution without any plot nor error test.
saddle_test_6.m: Error test of saddle_sparse.m for an individual element.
saddle_test_n.m: Error test of saddle_sparse.m for a n dimension problem printing full matrices to compare values.
saddle_test_n_2d.m: 2D saddle_sparse.m ideal test including error printing and 2D plots. The ground is at the bottom of the screen. You can change the examples by changing the initial wind conditions on lines 40-48 values. I would normally change the directions. For instance, these are some examples:
v0f(s)=[-1,1,0,0]; -> horizontal constant right flow
v0f(s)=[1,-1,0,0]; -> horizontal constant left flow
v0f(s)=[0,0,-1,1]; -> vertical constant up flow
v0f(s)=[0,0,-1,1]; -> vertical constant down flow
v0f(s)=[-1,0,1,0]; -> diagonal constant right-down flow
...
In the code commented, there is also an example with a hole.
saddle_test_n_3d.m: 3D saddle_sparse.m ideal test including error printing and 3D plots. The ideal case is a hill in the middle of the domain and a constant horizontal flow going into the mountain. You can change the terrain with the arguments of the add_terrain_to_mesh.m function on line 20. Also, you can change the initial flow in a very similar way to the saddle_test_n_2d.m on file sparse_assembly.m and line 51 (v0).
test_A.m: Runs the code for a similar 3D ideal case than saddle_test_n_3d.m but changing the number of elements. It interpolates to the same points using different resolutions. It takes a long time and a big memory to run. But, you can decrease the number of levels just to see what it does.
test_refinement.m: Tests the consistency of the results from test_A.m resulting Matlab file using error and 3D plots. So, you need to first run test_A.m. After that, you need to load the Matlab file created by test_A.m and run test_refinement.m.
test_WRF.m: Tests mass-consistent results initialized from an ideal case on WRF compared to the same WRF simulation after some timesteps printing errors, and 2d and 3D plots. I generated WRF_20.mat, WRF_45.mat, and WRF_241.mat to have light files with the important information to initialize the mass-consistent and run the comparisons. These three files are attached to this email and they correspond to 3 different examples that you can choose in line 5 of the test_WRF-m file. The three examples use the same initial wind flow (u0, v0, and w0), but they take different time steps to compare the mass-consistent to:
WRF_20.mat: timestep 20.
WRF_45.mat: timestep 45.
WRF_241.mat: timestep 241.
The other part of the code (the parent folder) solves the mass-consistent problem using multigrid instead of the saddlepoint problem and I did an excel file summarizing the different functions and explaining what they do. The file is on Onedrive under NASA/Angel/Projects/Massconsistent/inventory_quicwind.xlsx. Unfortunately, I didn't do something similar for the saddlepoint solution... However, in the same Onedrive folder, you can find also notes and some cases in png and fig format of some of the tests explained above. 

I need also to note that you need to first run the startup.m script from the parent folder before running the code since it uses some functions like big.


