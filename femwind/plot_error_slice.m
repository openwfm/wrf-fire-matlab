function plot_error_slice(e,r,X,tstring,params)
% plot_error_slice(e,r,F,K,X,tstring,params)
n = size(X{1});
nn=prod(n);
s=round(0.5+params.slice*n(2)*(1-eps));  % y index of slice
xx=squeeze(X{1}(:,s,:));
yy=squeeze(X{2}(:,s,:));
zz=squeeze(X{3}(:,s,:));
t=sprintf('slice %g y=%g %s',s,yy(1),tstring);
if ~isempty(r),
    residual=zeros(n);
    for i=1:nn
        [s1,s2,s3]=ind2sub(n,i);
        residual(s1,s2,s3)=r(i);
    end
    figure(params.res_slice_fig)
    mesh(xx,zz,squeeze(residual(:,s,:)));
    xlabel('horizontal')
    ylabel('vertical')
    ff = 'no_error_title';
    if isfield(params,ff) && ~getfield(params,ff)
        title(t)
    else
        title(['residual ',t])
    end
end
if ~isempty(e)
    lambda_err=zeros(n);
    for i=1:nn
        [s1,s2,s3]=ind2sub(n,i);
        lambda_err(s1,s2,s3)=e(i);
    end
    figure(params.err_slice_fig),clf
    l=squeeze(lambda_err(:,s,:));
    mesh(xx,zz,l)
    xlabel('horizontal')
    ylabel('vertical')
    title(['error ',t])
end
drawnow
