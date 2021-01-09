function [CX,CH]=center_mesh(X)
% CX=center_mesh(X)
% compute centers of mesh cells by averaging the coordinates of the cell corners
% in: 
%   X    {x,y,z}  3D arrays of coordinates of mesh corners
% out:
%   CX   {x,y,z}  3D arrays of coordinates of mesh cell centers
%   CW   {x,y,z}  3D arrays of coordinates of mesh horizontal face centers
%
CX = cell(1,3);
for k = 1:3
    % centers
    CX{k} = (X{k}(1:end-1,1:end-1,1:end-1)+X{k}(2:end,1:end-1,1:end-1)+...
            X{k}(1:end-1,2:end,1:end-1)+X{k}(1:end-1,1:end-1,2:end)+...
            X{k}(2:end,2:end,1:end-1)+X{k}(2:end,1:end-1,2:end)+...
            X{k}(1:end-1,2:end,2:end)+X{k}(2:end,2:end,2:end))/8;
end
% bottom
CB = X{3}(1:end-1,1:end-1,1)+X{3}(2:end,1:end-1,1)+...
     X{3}(1:end-1,2:end,1)+X{3}(1:end-1,1:end-1,1);
% center height above terrain
CH = zeros(size(CX{3}));
for k=1:size(CX{3},3)
     CH(:,:,k)=CX{3}(:,:,k)-CB;
end
end