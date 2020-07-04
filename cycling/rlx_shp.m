function out_shp = rlx_shp(in_shp,alpha,patch_size)
%does some kind of relaxation scheme on shapes
% local average over sub-matrix determined by patch_size
% alpha \in [0,1] weight of average

a = in_shp;
beta = 1-alpha;
new_shp = in_shp;
[n,m] = size(in_shp);
%patch_size = 4;
for i = patch_size+1:n-patch_size
    for j = patch_size+1:m-patch_size
        %a(i,j) = alpha*a(i,j)+beta/4*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1));
        patch = a(i-patch_size:i+patch_size,j-patch_size:j+patch_size);
        a(i,j) = alpha*a(i,j)+beta*mean(patch(:));
    end
end
out_shp = a;
end
