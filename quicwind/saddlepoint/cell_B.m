function B = cell_B(X)
% B=matrix_B(X)
% in
%   X     coordinates of the eight vertices of a cell (24 numbers)
% out
%   B     multiplication by B transforms fluxes through the faces of the
%   cell to 2 wind vectors in Cartesian coordinates
%   Xv    coordinates of two points where the wind vectors are defined
%   (later)

% v = B*u 
% u = flux through faces left x1, right x1, front x2, back x2, bottom x3,
% top x3 - always low to high coordinates
% v = x1,x2,x3,x1,x2,x3 cartesian components of the wind vector
% X{1} = x1 components of node coordinates, size 2x2x2, low-high x1,x2,x3
% X{2} = x2 components of node coordinates, size 2x2x2, low-high x1,x2,x3
% X{3} = x3 components of node coordinates, size 2x2x2, low-high x1,x2,x3
figure, plot_mesh_3d(X);

end