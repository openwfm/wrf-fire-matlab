function F=ndt_f_assembly_fortran(A,X,u0, lambda, params)
% call fortran version

[~,F,~]= sparse_assembly(A,X,u0,lambda,params);

%if params.test_fortran
if 1==1
    disp('testing if sprase assembly of F same result in fortran')
    exe = './fortran/ndt_f_test.exe';
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
        
        F_fort=swap23(read_array_nd('Fvec'));
        err= norm(F(:)-F_fort(:),inf)
        
        if abs(err)>1e-6
            F_fort;
            error('sprase assembly of F different')
        end
        disp('sprase assembly of F result in fortran OK')
    else
        warning(['file ',exe,' does not exist'])
    end
end
end
