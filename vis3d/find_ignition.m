function r=find_ignition(f)
% r=find_ignition(f)
% finds the first frame with nodes where fire_area > 0
% arguments:
%     f wrfout file name
% returns:
%     structure with ignited nodes and their coordinates

p = nc2struct(f,{'Times'},{});
times=char(p.times)';
steps=size(times,1);
r = [];
for step=1:steps
    p=nc2struct(f,{'FIRE_AREA'},{},step)
    [i,j,a]=find(p.fire_area);
    n = length(i);
    disp([num2str(step),' ',times(step,:),' ignited ',num2str(n)])
    if n
        q=nc2struct(f,{'FXLONG','FXLAT'},{},step);
        r.file=f;
        r.frame=step;
        r.i=i';
        r.j=j';
        r.times=times(step,:);
        for k=n:-1:1
            r.fxlong(k)=q.fxlong(i(k),j(k));            
            r.fxlat(k)=q.fxlat(i(k),j(k));            
            r.fire_area(k)=p.fire_area(i(k),j(k));
        end
        break
    end
end
end
        
        
        


