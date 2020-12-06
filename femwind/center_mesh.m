function CX=center_mesh(X)
% CX=center_mesh(X)
% compute centers of mesh cells by averaging the coordinates of the cell corners
% in: 
%   X    {x,y,z}  3D arrays of coordinates of mesh corners
% out:
%   DX   {x,y,z}  3D arrays of coordinates of mesh cell centers 
CX = cell(1,3);
for k = 1:3
    CX{k} = (X{k}(1:end-1,1:end-1,1:end-1)+X{k}(2:end,1:end-1,1:end-1)+X{k}(1:end-1,2:end,1:end-1)+X{k}(1:end-1,1:end-1,2:end)+X{k}(2:end,2:end,1:end-1)+X{k}(2:end,1:end-1,2:end)+X{k}(1:end-1,2:end,2:end)+X{k}(2:end,2:end,2:end))/8;
end
end