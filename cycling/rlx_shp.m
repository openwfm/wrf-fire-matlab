function out_shp = rlx_shp(in_shp,alpha,mask)
%does some kind of relaxation scheme on shapes

a = in_shp;
beta = 1-alpha;
new_shp = in_shp;
[n,m] = size(in_shp);
for i = 2:n-1
    for j = 2:m-1
        a(i,j) = alpha*a(i,j)+mask(i,j)*beta/4*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1));
    end
end
out_shp = a;
end
