function [h,u,v,names]=prof(p,long,lat)
names={'layer1','layer2','layer3','layer4','can_top','fwh'};
m=4;
for i=1:m
    h(i)=interps(p.xlong,p.xlat,p.height(:,:,i),long,lat); 
    u(i)=interps(p.xlong,p.xlat,p.uc(:,:,i),long,lat); 
    v(i)=interps(p.xlong,p.xlat,p.vc(:,:,i),long,lat); 
end
i=m+1;
h(i)=interps(p.fxlong,p.fxlat,p.can_top,long,lat); 
u(i)=interps(p.fxlong,p.fxlat,p.cuf,long,lat); 
v(i)=interps(p.fxlong,p.fxlat,p.cvf,long,lat);
i=m+2;
h(i)=interps(p.fxlong,p.fxlat,p.fwh,long,lat); 
u(i)=interps(p.fxlong,p.fxlat,p.uf,long,lat); 
v(i)=interps(p.fxlong,p.fxlat,p.vf,long,lat); 

names,h,u,v
end

function vq=interps(x,y,v,xq,yq)
   F = scatteredInterpolant(x(:),y(:),v(:));
   vq=F(xq,yq);
end
 

