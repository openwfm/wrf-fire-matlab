function err=wind2flux_trans_test
nx=50; ny=30; nz=10;
h=rand(1,3);
X = regular_mesh([nx,ny,nz],h,1.2);
X = add_terrain_to_mesh(X,'hill','squash',0.1);

% random test winds
U = grad3z(rand(size(X{1})-1),[1 1 1]); 
V = grad3z(rand(size(X{1})-1),[1 1 1]);

% Test that <Mu,v>=<u,M^Tv>
rhs = aprod3(wind2flux(U,X),V);
lhs = aprod3(U,wind2flux_trans(V,X));
err = big(rhs-lhs);
end

function a=aprod(x,y)
    a = dot(x(:),y(:));
end

function a=aprod3(x,y)
    a = aprod(x{1},y{1})+aprod(x{2},y{2})+aprod(x{3},y{3});
end