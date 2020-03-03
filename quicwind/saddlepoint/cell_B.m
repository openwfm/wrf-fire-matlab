function B = cell_B(dx,dy,z,i,j,k)
% B=matrix_B(X)
% in
%   dx    spacing of the grid in the x direction
%   dy    spacing of the grid in the y direction
%   z     matrix of the height of the nodes
% out
%   B     multiplication by B transforms fluxes through the faces of the
%   cell to 2 wind vectors in Cartesian coordinates
%   Xv    coordinates of two points where the wind vectors are defined
%   (later)

% v = B*u 
% u = flux through faces left x1, right x1, front x2, back x2, bottom x3, top x3 - always low to high coordinates
% v = v1(1),v2(1),v1(2),v2(2),v1(3),v2(3) cartesian components of the wind vectors
% X{1} = x1 components of node coordinates, size 2x2x2, low-high x1,x2,x3
% X{2} = x2 components of node coordinates, size 2x2x2, low-high x1,x2,x3
% X{3} = x3 components of node coordinates, size 2x2x2, low-high x1,x2,x3

%figure, plot_mesh_3d(X);
% tangents
txk = (z(i+1,j,k)-z(i,j,k)+z(i+1,j+1,k)-z(i,j+1,k))/(2*dx); % tangent of thetax at k
tyk = (z(i,j+1,k)-z(i,j,k)+z(i+1,j+1,k)-z(i+1,j,k))/(2*dy); % tangent of thetay at k
txk1 = (z(i+1,j,k+1)-z(i,j,k+1)+z(i+1,j+1,k+1)-z(i,j+1,k+1))/(2*dx); % tangent of thetax at k+1
tyk1 = (z(i,j+1,k+1)-z(i,j,k+1)+z(i+1,j+1,k+1)-z(i+1,j,k+1))/(2*dy); % tangent of thetay at k+1
% secants 
sck = sqrt((1+txk*txk)*(1+tyk*tyk)); % sec(thetax)*sec(thetay) at k 
sck1 = sqrt((1+txk1*txk1)*(1+tyk1*tyk1)); % sec(thetax)*sec(thetay) at k+1
% areas at the faces
a1 = dy*.5*(z(i,j,k+1)-z(i,j,k)+z(i,j+1,k+1)-z(i,j+1,k)); % left face area
a2 = dy*.5*(z(i+1,j,k+1)-z(i+1,j,k)+z(i+1,j+1,k+1)-z(i+1,j+1,k)); % right face area
a3 = dx*.5*(z(i,j,k+1)-z(i,j,k)+z(i+1,j,k+1)-z(i+1,j,k)); % front face area 
a4 = dx*.5*(z(i,j+1,k+1)-z(i,j+1,k)+z(i+1,j+1,k+1)-z(i+1,j+1,k)); % front face area 
a5 = (dx*dy)*sck; % bottom face area 
a6 = (dx*dy)*sck1; % top face area
a = [a1,a2,a3,a4,a5,a6];
Binv = diag(a)* ...
        [-1,0,0,0,0,0;       % 1: flux through the left face from v1
         0,1,0,0,0,0;        % 2: flux through the right face from v2
         0,0,-1,0,0,0;       % 3: flux through the front face from v1
         0,0,0,1,0,0;        % 4: flux through the back face from v2
         txk,0,tyk,0,-1,0;   % 5: flux through the bottom face from v1
         0,-txk1,0,-tyk1,0,1]; % 6: flux through the top face from v2
%disp('columns are: left, right, front, back, bottom, and top');
%disp('rows are: v1(1), v2(1), v1(2), v2(2), v1(3), and v2(3)');
B = inv(Binv);
end