% a quick script to prepare terrain height for fewmind 
% from elevation in a geotiff

filename='GranitePeak-DEM.tif';
sr_x=10
sr_y=10
disp(['reading ',filename])
t = Tiff(filename,'r');
info = geotiffinfo(filename)
dx = info.PixelScale(1)
dy = info.PixelScale(2)
d = read(t);
ds = size(d);
disp(['grid size ',num2str(size(d)),' min ',num2str(min(d(:))),' max ',num2str(max(d(:)))])
isza = floor(ds(1)/sr_x)
iszf = isza * sr_x
jsza = floor(ds(2)/sr_y)
jszf = jsza * sr_y
zsf = d(1:iszf,1:jszf);
[xf,yf] = ndgrid(1:iszf,1:jszf);
[xa,ya] = ndgrid(([1:isza]-0.5)*sr_x,([1:jsza]-0.5)*sr_y);

ht1 = interp2(1:jszf,1:iszf,zsf,ya(:),xa(:));
if isnan(sum(ht1)),error('NaN in interpolant'),end
ht = reshape(ht1,size(xa));

figure(1), mesh(xf,yf,zsf)
%figure(2), 
hold on
mesh(xa,ya,ht)
hold off 
drawnow


fzsf='input_zsf';
disp(['writing ',fzsf])
write_array_2d(fzsf,zsf)

fht='input_ht';
disp(['writing ',fht])
write_array_2d(fht,ht)

