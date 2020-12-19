disp('femwind_test')

do_plot=0

scale=[1,2,4,8,16]
% scale=[1]
for icase = 1:length(scale)
sc = scale(icase)    

% dimensions in elements
% sc=1; % mesh refinement scale
n = sc*[20,10,5];
sc2=1;
n(1:2)=n(1:2)*sc2
h = [10,10,10]/sc;
fprintf('linear array of %ix%ix%i cells\n',n(1),n(2),n(3))
da=[1 1 1]
string_diag_A=sprintf('%g %g %g',da);
A = diag(da);
lambda = zeros(prod(n+1),1); % placeholder solution

% creating the grid
expand=1.0
X = regular_mesh(n,h,expand^(1/sc));
X = add_terrain_to_mesh(X,'hill','squash',0.3);
CX = center_mesh(X);

% initial wind at the centers of the elements
U0={ones(n),zeros(n),zeros(n)};

if 1

    % show mesh
    hold off
    figure(1)
    plot_mesh_3d(X), hold on, 
    title('The wind mesh, wind vector in centers, lambda in corners')

    % show initial wind
    figure(2)
    plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
    plot_wind_3d(CX,U0)
    hold off
    axis equal
    title('Initial wind')

    % show initial wind
    figure(3)
    plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
    plot_wind_3d(CX,U0,1)
    hold off
    axis equal
    title('Initial wind lowest layer')
end

% assemble sparse system matrix
[K,F,~] = sparse_assembly(A,X,U0,lambda);

% dirichlet boundary conditions
[K,F]=apply_boundary_conditions(K,F,X);

% solve the equations
% [lambda,it] = sparse_solve(K,F,X,'s');
[lambda,it] = sparse_solve(K,F,X,'2');

% assemble final wind
[~,~,W] = sparse_assembly(A,X,U0,lambda);

if do_plot
    % plot resulting wind
    figure(4)
    plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
    plot_wind_3d(CX,W)
    hold off
    axis equal
    title(['Final wind a=',string_diag_A])

    figure(5)
    plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
    plot_wind_3d(CX,W,1)
    hold off
    axis equal
    title(['Final wind lowest layer a=',string_diag_A])

    figure(6)
    plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
    plot_wind_3d(CX,W,1:2)
    hold off
    axis equal
    title(['Final wind lowest layers a=',string_diag_A])
    
end


    figure(7)
    height=10;
    [XH,WH]=wind_at_h(X,CX,W,[20,20,1],...
        [min(X{1}(:)),max(X{1}(:)),...
        min(X{3}(:)),max(X{3}(:)),...
        height,height]); 
    plot_wind_3d(XH,WH)
    hold on
    plot_mesh_3d(X,[1,n(1),1,n(2),1,1])
    hold off
    axis equal
    title(['Final wind with a=',string_diag_A,' at ',num2str(height),' above terrain'])

    
    s(icase).sc=sc;
    s(icase).n=n;
    s(icase).X=X;
    s(icase).CX=CX;
    s(icase).XH=XH;
    s(icase).h=h;
    s(icase).WH=WH;
    s(icase).W=W;
    s(icase).A=A;
    s(icase).U0=U0;
    s(icase).it=it;
    save s.mat s
    
end
% condition_number=scond(K)
n
