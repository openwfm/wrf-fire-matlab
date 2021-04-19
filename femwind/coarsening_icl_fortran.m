function [hzc,icl_z]=coarsening_icl_fortran(dx,dy,dz,params)
% call fortran version

[hzc,icl_z]=coarsening_icl(dx,dy,dz,params)

if params.test_fortran
    disp('testing if coarsening_icl_test same result in fortran')
    exe = './fortran/coarsening_icl_test.exe';
    if exist(exe,'file') 
        %Writing all arrays to text files for use by fortran tester
        write_array(dx,'dx');
        write_array(dy,'dy');
        write_array(dz,'dz');
        write_array(diag(params.a),'A');
        write_array(params.minaspect,'minaspect');
        write_array(params.maxaspect,'maxaspect');
        system(exe)
        hzcf=[read_array('cr_x'),read_array('cr_y')];
        icl_zf = read_array('icl_z'); icl_zf=icl_zf(:)'; % row
        if any(hzc ~= hzcf) || any(icl_zf ~= icl_z)
            hzc,icl_z
            hzcf,icl_zf
            error('coarsening_icl_fortran different')
        end
        disp('coarsening_icl_test same result in fortran OK')
    else
        warning(['file ',exe,' does not exist'])
    end
end
end
