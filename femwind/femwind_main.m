function params=femwind_main(params)
% usage:
% params=femwind_main    output params structure and exit
% femwind_test(params)   run with params specified and copy to output
% 
%           *** WARNING ***
% DO NOT EDIT THIS FILE to set params or add outputs.
% Instead, call the function to get the params, add or delete any fields you need, and call it again.
% Debugging outputs can be added to the params structure as additional field.
%
%   *** all_test must pass fter any edits and before any merge ***
%
%           *** WARNING ***


if ~exist('params','var') | isempty(params)
    params = params_defaults
    return 
end

if params.save_files >= 0
    diary(['femwind_',params.save_file_prefix,'_diary.txt'])
end
disp('femwind_main')
format compact

if isfield(params,'mesh_top')
    if params.mesh_top>0
        disp('given params.mesh_top>0, computing params.expand') 
        nz = params.nelem3(3); % elements in the vertical direction
        a = params.mesh_top/params.h(3); % desired height as multiple of first layer
        params.expand = findq(a,nz);
    end

end

params


for sc2 = params.sc2_all
    for sc = params.sc_all
        sc,sc2
        
        disp('setting up test case')

        nel = sc*params.nelem3;  % elements in the 3 directions
        nel(1:2)=nel(1:2)*sc2
        h = params.h/sc;
        fprintf('mesh of %ix%ix%i cells\n',nel(1),nel(2),nel(3))
        params.id=sprintf('%ix%ix%i',nel); % to pass around 
        string_diag_A=sprintf('%g %g %g',params.a); % for figure titles
        A = diag(params.a.^2);

        % creating the grid
        expand=params.expand;
        X = regular_mesh(nel,h,params.expand^(1/sc));
        X = add_terrain_to_mesh(X,...
            params.terrain_shape,params.terrain_top,params.terrain_height);
        [CX,CH] = center_mesh(X); % get midpoints of elements

        % initial wind at the centers of the elements
        rng(1);
        switch params.initial_wind
            case 'uniform'
                disp('initial wind uniform in x direction')
                U0={ones(nel),zeros(nel),zeros(nel)};
            case 'random-z'
                % to test iterative methods with non-smooth initial error
                disp('initial wind uniform in x direction random in z direction')
                U0={ones(nel),zeros(nel),randn(nel)};
            case 'random-xz'
                % to test iterative methods with non-smooth initial error
                disp('initial wind uniform in x direction and z direction')
                U0={ones(nel),zeros(nel),ones(nel)};
            case 'log'
                disp('initial wind log profile in x direction')
                U0={log(max(1,CH/params.roughness_height)),zeros(nel),zeros(nel)};
            otherwise
                error(['unknown initial wind ',params.initial_wind])
        end
        if params.graphics>0
            disp('graphics: problem setup')
            % show mesh
            figure(1),clf
            plot_mesh_3d(X)
            axis equal
            title('The wind mesh, wind vector in centers, lambda in corners')
        end

        if params.graphics>0
            % show initial wind
            figure(2),clf
            plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,1,1]), hold on, 
            plot_wind_3d(CX,U0)
            hold off
            axis equal
            title('Initial wind')
        end

        if params.graphics>1
            % show initial wind
            figure(3),clf
            plot_mesh_3d(X,[1,nel(1),1,nel(2)+1,1,1]), hold on, 
            plot_wind_3d(CX,U0,1)
            hold off
            axis equal
            title('Initial wind lowest layer')
        end
        
%         Plot initial streamlines
        if params.graphics > 2
            figure(4), clf
            plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,2,2])
            hold on
            wind_streamlines(X, CX, U0, params);
            hold off
        end
        diary; diary

        % solve
        if params.run_fortran && params.run_matlab
                [W,rate(sc,sc2)]=femwind_fortran(A,X,U0,params);
                [Wm,rate(sc,sc2)]=femwind_solve(A,X,U0,params);
                for i=1:3, err_compare(i)=big(W{i}-Wm{i})/big(W{i});end
                fprintf('femwind_fortran and femwind_solve done\n')
                err=max(err_compare);
                fprintf('relative max error matlab vs fortran %g\n',err)
                if err>1e-4,
                    warning('error too large')
                end
        elseif params.run_fortran
                [W,rate(sc,sc2)]=femwind_fortran(A,X,U0,params);
                fprintf('femwind_fortran done ')
        else
                [W,rate(sc,sc2)]=femwind_solve(A,X,U0,params);
                fprintf('femwind_solve done ')
        end
        fprintf('rate(%i,%i)=%g\n',sc,sc2,rate(sc,sc2))
        
        if params.graphics>1
            disp('graphics: solution')

            % plot resulting wind
            figure(5),crate(sc,sc2)
            plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,1,1]), hold on, 
            plot_wind_3d(CX,W)
            hold off
            axis equal
            title(['Final wind a=',string_diag_A])

            figure(6),clf
            plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,1,1]), hold on, 
            plot_wind_3d(CX,W,1)
            hold off
            axis equal
            title(['Final wind lowest layer a=',string_diag_A])

            figure(7),clf
            plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,1,1]), hold on, 
            plot_wind_3d(CX,W,1:2)
            hold off
            axis equal
            title(['Final wind lowest layers a=',string_diag_A])

        end

        if params.graphics>0
            disp('graphics: wind_at_h')

            figure(8),clf
            height=10;
            bbox=[min(X{1}(:)),max(X{1}(:)),...
                min(X{2}(:)),max(X{2}(:)),...
                height,height];
            [XH,WH]=wind_at_h(X,CX,W,[20,20,1],bbox);
            plot_wind_3d(XH,WH)
            hold on
            plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,1,1])
            hold off
            axis equal
            title(['Final wind with a=',string_diag_A,' at ',num2str(height),' above terrain'])
        end
        if params.graphics > 2
            figure(9),clf
            plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,2,2])
            hold on
            wind_streamlines(X, CX, W, params)
            hold off
        end    
        params.rate=rate;
        if params.save_files>1,
            sfile=[params.save_file_prefix,'_workspace.mat'];
            disp(['saving workspace to ',sfile])
            save(sfile,'-v7.3') 
        end
        if params.save_files>0,
            sfile=[params.save_file_prefix,'_params.mat'];
            disp(['saving params to ',sfile])
            save(sfile,'params','-v7.3') 
        end
    end
end
