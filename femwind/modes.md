commit 4307bbffee2e537e9ed71922b849ab3586cf3769
Author: Jan Mandel <jan.mandel@gmail.com>
Date:   Sun Dec 27 14:44:29 2020 -0700

params 
         graphics: 1
            expand: 1.2000
                sc: [1 2]
               sc2: [1 2 4 8 16 32]
            nelem3: [20 20 8]
                 h: [10 10 10]
                da: [1 1 1]
      initial_wind: 1
     terrain_shape: 'hill'
       terrain_top: 'squash'
    terrain_height: 0.2000
            solver: '2-level'
             maxit: 50
        coarsening: '2 linear'
         smoothing: 'vertical lines'
           nsmooth: 5
            restol: 1.0000e-06
             exact: 0
             slice: 0.5000
     err_slice_fig: 12
     res_slice_fig: 13
         maxaspect: 3

Rate
 Columns 1 through 7

    0.0605    0.0947         0    0.0925         0         0         0
    0.0494    0.0804         0    0.0844         0         0         0

  Columns 8 through 14

    0.1028         0         0         0         0         0         0
    0.1110         0         0         0         0         0         0

  Columns 15 through 21

         0    0.1085         0         0         0         0         0
         0    0.1152         0         0         0         0         0

  Columns 22 through 28

         0         0         0         0         0         0         0
         0         0         0         0         0         0         0

  Columns 29 through 32

         0         0         0    0.1104
         0         0         0         0

**********************************************************************

commit 660e41ab77cfc94714c6056b9af916dc0f878a07 (HEAD -> femwind)
Author: Jan Mandel <jan.mandel@gmail.com>
Date:   Wed Dec 30 19:52:24 2020 -0700

params = 
          graphics: 2
            expand: 1.2000
                sc: 1
               sc2: 1
            nelem3: [21 21 8]
                 h: [10 10 10]
                da: [1 1 1]
      initial_wind: 1
     terrain_shape: 'hill'
       terrain_top: 'squash'
    terrain_height: 0.2000
            solver: '2-level'
             maxit: 50
        coarsening: '2 linear'
         smoothing: 'vertical lines'
           nsmooth: 3
            restol: 1.0000e-06
             exact: 0
             slice: 0.5000
     err_slice_fig: 12
     res_slice_fig: 13
    iterations_fig: 14
         maxaspect: 3
            levels: 3
    nsmooth_coarse: 2
      maxit_coarse: 8

rate =
    0.1335    0.2641

same, with params.levels=2:
rate =
    0.1335    0.1763

************************************************************
commit 3306348da63093b8e4500cb7960d9af5c1c836f7
Author: Jan Mandel <jan.mandel@gmail.com>
Date:   Tue Dec 29 22:54:34 2020 -0700

params = 

  struct with fields:

          graphics: 1
            expand: 1.2000
                sc: [1 2]
               sc2: [1 2 4 8 16 32]
            nelem3: [22 22 8]
                 h: [10 10 10]
                da: [1 1 1]
      initial_wind: 1
     terrain_shape: 'hill'
       terrain_top: 'squash'
    terrain_height: 0.2000
            solver: '2-level'
             maxit: 50
        coarsening: '2 linear'
         smoothing: 'vertical sweeps'
           nsmooth: 3
            restol: 1.0000e-06
             exact: 0
             slice: 0.5000
     err_slice_fig: 12
     res_slice_fig: 13
    iterations_fig: 14
         maxaspect: 3
            levels: 3
    nsmooth_coarse: 2
      maxit_coarse: 8

  Columns 1 through 13

    0.1937    0.2442         0    0.4631         0         0         0    0.5760         0         0         0         0         0
    0.2156    0.2512         0    0.2892         0         0         0    0.3404         0         0         0         0         0

  Columns 14 through 26

         0         0    0.6109         0         0         0         0         0         0         0         0         0         0
         0         0         0         0         0         0         0         0         0         0         0         0         0

  Columns 27 through 32

         0         0         0         0         0    0.6228
         0         0         0         0         0         0

reason: slow coarse convergence

after fixing coarse:

commit ba01f53f1e4bddce83f2e90b18e477eb0443c4d6
Author: Jan Mandel <jan.mandel@gmail.com>
Date:   Thu Dec 31 12:39:07 2020 -0700

params = 
          graphics: 1
            expand: 1.2000
            sc_all: [1 2]
           sc2_all: [1 2 4 8 16 32]
            nelem3: [22 22 8]
                 h: [10 10 10]
                da: [1 1 1]
      initial_wind: 1
     terrain_shape: 'hill'
       terrain_top: 'squash'
    terrain_height: 0.2000
            solver: '2-level'
             maxit: 50
        coarsening: '2 linear'
         smoothing: 'vertical sweeps'
           nsmooth: 3
            restol: 1.0000e-06
             exact: 0
             slice: 0.5000
     err_slice_fig: 12
     res_slice_fig: 13
    iterations_fig: 14
         maxaspect: 3
            levels: 3
    nsmooth_coarse: 2
      maxit_coarse: 8


rate =

  Columns 1 through 11

    0.1919    0.2300         0    0.2309         0         0         0    0.2629         0         0         0
    0.2148    0.2436         0    0.2516         0         0         0         0         0         0         0

  Columns 12 through 22

         0         0         0         0    0.2665         0         0         0         0         0         0
         0         0         0         0         0         0         0         0         0         0         0

  Columns 23 through 32

         0         0         0         0         0         0         0         0         0    0.2667
         0         0         0         0         0         0         0         0         0         0
