function [h,u,v,names]=prof(p,long,lat)
m=4;
long
lat

disp('nearest in xlong,xlat')
[i,j]=minij(abs(p.xlong-long)+abs(p.xlat-lat))
lo=p.xlong(i,j)
la=p.xlat(i,j)

for k=1:m
    k,
    layer_height=p.height(i-1:i+1,j-1:j+1,k)
    layer_uc=p.uc(i-1:i+1,j-1:j+1,k)
    layer_vc=p.vc(i-1:i+1,j-1:j+1,k)
end

disp('nearest in fxlong,fxlat')
[fi,fj]=minij(abs(p.fxlong-long)+abs(p.fxlat-lat))
flo=p.fxlong(i,j)
fla=p.fxlat(i,j)
cuf=p.cuf(fi,fj)
cvf=p.cvf(fi,fj)



names={'layer1','layer2','layer3','layer4','can_top','uf/vf','uah/vah'};
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
i=m+3;
% h(i)=interps(p.xlong,p.xlat,p.fwh,long,lat); 
ua=(p.uah(1:end-1,:)+p.uah(2:end,:))/2;
va=(p.vah(:,1:end-1)+p.vah(:,2:end))/2;
u(i)=interps(p.xlong,p.xlat,ua,long,lat); 
v(i)=interps(p.xlong,p.xlat,va,long,lat); 

names,h,u,v
end

function vq=interps(x,y,v,xq,yq)
   F = scatteredInterpolant(x(:),y(:),v(:));
   vq=F(xq,yq);
end
 
function [i,j]=minij(a)
   [ix,jx,v]=find(a);
   [~,k]=min(v);
   i=ix(k);
   j=jx(k);
end
  

