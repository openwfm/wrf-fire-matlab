function [w,rate]=femwind_solve_fortran(A,X,u0,params)

exe  = './fortran/femwind_solve_test.exe';
if exist(exe,'file') & params.run_fortran
    write_array_nd(swap23(X{1}),'X');
    write_array_nd(swap23(X{2}),'Y');
    write_array_nd(swap23(X{3}),'Z');
    write_array_nd(swap23(u0{1}),'Xu0');
    write_array_nd(swap23(u0{2}),'Yu0');
    write_array_nd(swap23(u0{3}),'Zu0');

    write_array(A,'A')
    write_array(iflags,'iflags')

    % defaults
    Xw =zeros(nel);
    Yw =zeros(nel);
    Zw =zeros(nel);
    rate=0;

    try
        disp(['running ',exe])
        system(exe);
        Xw=swap23(read_array_nd('Xw'));
        Yw=swap23(read_array_nd('Yw'));
        Zw=swap23(read_array_nd('Zw'));
        rate=read_array_nd('rate');
    catch
        disp([exe,' failed'])
    end

    wf = {Xw,Yw,Zw};
    w = wf;
end

if params.run_matlab
    
    [wm,rate]=femwind_solve(A,X,u0,params)
    w = wm;
end

if exist('wf','var') & exist('wm','var')

    tol = 10*eps(single(1.));
    
    wwf=cell2mat(uf);
    wwm=cell2mat(um);
 
    err= norm(wwf(:)-wwm(:),inf)/norm(wwm(:));
    if err < tol
        fprintf('error %g OK, tol = %g\n',err,tol)
    else
        warning(sprintf('error %g too large, tol=%g',err,tol))
    end

end
