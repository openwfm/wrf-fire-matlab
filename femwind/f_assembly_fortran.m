function F=f_assembly_fortran(A,X,u0, lambda, params)
% call fortran version

F= f_assembly(A,X,u0);

%if params.test_fortran
if 1==1
    disp('testing if sparse assembly of F same result in fortran')
    exe = './fortran/f_assembly_test.exe';
    if exist(exe,'file') 
        %Writing all arrays to text files for use by fortran tester
        write_array_nd(swap23(X{1}),'X');
        write_array_nd(swap23(X{2}),'Y');   
        write_array_nd(swap23(X{3}),'Z');
        write_array_nd(swap23(u0{1}),'Xu0');
        write_array_nd(swap23(u0{2}),'Yu0');
        write_array_nd(swap23(u0{3}),'Zu0');
        write_array(A,'A');
        
        system(exe)
        
        F_fort=swap23(read_array_nd('F'));
        err= norm(F(:)-F_fort(:),inf)
        
        tol = 10*eps(single(1.));
        if err < tol
            fprintf('error %g OK, tol = %g\n',err,tol)
        elseif err < 0.2
            warning(sprintf('error %g too large\n',err))
        else
            error(sprintf('error %g too large, tol=%g',err,tol))
        end

    else
        warning(['file ',exe,' does not exist'])
    end
end
end
