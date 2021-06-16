function plot_error_slice(e,r,X,tstring,params)
% plot_error_slice(e,r,F,K,X,tstring,params)
if params.graphics < 0
    return
end
n = size(X{1});
nn=prod(n);
s=round(0.5+params.slice*n(2)*(1-eps));  % y index of slice
xx=squeeze(X{1}(:,s,:));
yy=squeeze(X{2}(:,s,:));
zz=squeeze(X{3}(:,s,:));
t=sprintf('slice %g y=%g %s',s,yy(1),tstring);
ff = 'no_error_title';
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
    if isfield(params,ff) && ~getfield(params,ff)
        title(t)
    else
        title(['error ',t])
    end
    
end
if params.graphics > 2
    zoom_center = false;
    figure(100)
    residual=zeros(n);
    if zoom_center
        cx = (max(max(xx))-min(min(xx)))/2;
        mesh(xx,zz,squeeze(residual(:,s,:))); view(2); xlim([min(min(xx)),max(max(xx)),cx-500,cx+500]);
    else
        mesh(xx,zz,squeeze(residual(:,s,:))); view(2); xlim([min(min(xx)),max(max(xx))]);
    end
    nslices=size(dir('./slice_l*.png'),1);
    saveas(gcf,sprintf('slice_l%02d.png',nslices+1))
end
drawnow
end