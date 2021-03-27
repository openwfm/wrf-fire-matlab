function check_P(P,X,XC)
n = size(X{1});
nn=prod(n);
nc=size(XC{1});
YC = XC{3};
for i=nc(3):-1:1
    YC(:,:,i)=YC(:,:,i)-YC(:,:,1);
end
for i=1:3; 
    xmax(i)=max(XC{i}(:));
    xmin(i)=min(XC{i}(:));
end
f = @(x1,x2,y)(y-xmin(3))/(xmax(3)-xmin(3));
FC = f(XC{1},XC{2},YC);
F = reshape(P*FC(:),n);
p.no_error_title=1;
p.slice=0.5;
p.err_slice_fig=1;
plot_error_slice(F,[],X,'interpolation',p) 
p.err_slice_fig=2;
plot_error_slice(FC,[],XC,'coarse',p) 
end