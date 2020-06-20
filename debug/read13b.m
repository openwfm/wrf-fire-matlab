vars={'UAH','VAH','U_2','V_2','F_LINEINT','F_LINEINT2','F_ROS','F_ROS13','F_ROS21','F_ROS23','F_ROS31','F_ROS32','F_ROS33','F_ROSX','F_ROSY','UF','VF'}
p10=nc2struct('wrf-10x10.2/wrfout_d03_2020-05-20_15:00:00',vars,{},2)
p9=nc2struct('wrf-9x9.2/wrfout_d03_2020-05-20_15:00:00',vars,{},2)
