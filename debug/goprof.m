 % debugging canopy interpolation
 step=100
 file='em_rxcadre2_devel/wrfout_d05_2012-11-11_12:00:00'
 long=-86.73
 lat=30.53
 p=wrfatm2struct(file,step);
 pp=nc2struct(file,{'CAN_TOP','CUF','CVF','FWH','FZ0'},{},step);
 p.can_top=pp.can_top;
 p.cuf=pp.cuf;
 p.cvf=pp.cvf;
 p.fwh=pp.fwh;
 p.fz0=pp.fz0;
 prof(p,long,lat);