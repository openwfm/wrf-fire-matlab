% file in: /home/math/farguella/wrf-fire-matlab/quicwind/saddlepoint/test_A.mat
% produced by running: test_A.m
load('test_A.mat')

n = [8,8,8];
h = [.5,.5,.5];
levels = 4;

X = regular_mesh(n,h,1);
thx = h(1)*[0:n(1)]'*ones(1,n(2)+1);
X = add_terrain_to_mesh(X,thx,'shift');
X = add_terrain_to_mesh(X,'hill','squash',0.1);
CX = cell(1,3);
for k = 1:3
    CX{k} = (X{k}(1:end-1,1:end-1,1:end-1)+X{k}(2:end,1:end-1,1:end-1)+X{k}(1:end-1,2:end,1:end-1)+X{k}(1:end-1,1:end-1,2:end)+X{k}(2:end,2:end,1:end-1)+X{k}(2:end,1:end-1,2:end)+X{k}(1:end-1,2:end,2:end)+X{k}(2:end,2:end,2:end))/8;
end
xx = CX{1}; yy = CX{2}; zz = CX{3};

for r=1:levels
    nr = [4,4,4]*2^(r-1);
    U{r} = F{r,1}(CX{1},CX{2},CX{3});
    V{r} = F{r,2}(CX{1},CX{2},CX{3});
    W{r} = F{r,3}(CX{1},CX{2},CX{3});
    figure(r), 
    plot_mesh_3d(X), 
    hold on, 
    quiver3(xx(:),yy(:),zz(:),U{r}(:),V{r}(:),W{r}(:),'LineWidth',2), 
    hold off,
    xlabel('x'), ylabel('y'), zlabel('z'), 
    title(['Mass-consistent wind ',num2str(nr)]);
end

%levels = 3;
for r=1:levels-1
    Udiff{r} = U{levels}-U{r};
    Vdiff{r} = V{levels}-V{r};
    Wdiff{r} = W{levels}-W{r};
    diff = [big(Udiff{r}),big(Vdiff{r}),big(Wdiff{r})],
    
    figure(r+4), 
    plot_mesh_3d(X), 
    hold on, 
    quiver3(xx(:),yy(:),zz(:),Udiff{r}(:),Vdiff{r}(:),0*Wdiff{r}(:),'LineWidth',2), 
    hold off,
    xlabel('x'), ylabel('y'), zlabel('z'),
    title(['Difference mass-consistent wind level ', num2str(levels),' minus level ',num2str(r),' with difference ',num2str(diff)]);
end