function [v,u] = rotate_cell(uIn)
% The function rotates and the flips over the off-diagaonal a 2x2 cell so that the
% largest value is in the lower left corner and largest valu on the
% diagonal is on the upper left corner.
% Usage:
% v = rotate_cell(uIn) will rotate permute the rows of matrix Uin.
% [v,u] = rotate_cell(~,n) will generate a nx4 matrix with random
% perumations of the row [1 2 3 4] and them permute them to correct order.
% Used for diagnostics.
% Inputs:
%    uIn - nx4 matrix whose rows will be permuted. Passing an integer as
%    uIn will create testing values with size uIn x 4
% Outputs:
%    v - matrix with rows to be permuted
%    u - optional matrix with created rows or original matrix.


% square grid with
%       t01 ----- t11
%        |         |
%        |         |
%       t00 ----- t10

% Optionally create matrix with rows [t00 t01 t10 t11]
% as perumatations of [4 3 2 1]
if length(uIn) == 1
    n = uIn;
    u = zeros(uIn,4);
    for i = 1:n
        u(i,:) = randperm(4);
    end
    v = u;
else
    v = uIn;
    u = uIn;
    [n,~] = size(v);
end

%rearrange so maximum value is in the first position, always.
%then put largest of diagonals in upper left.
for i = 1:n
    %find location of largest element
    [~,idx] = max(v(i,:));
    
    %temp row vector
    vt = v(i,:);
    %rotate to put largest first
    if idx ~= 1
        if idx == 2
            vt = vt([2 4 1 3]);
        elseif idx == 3
            vt = vt([3 1 4 2]);
        elseif idx == 4
            vt = vt([4 1 2 3]);
        end
        
    end
    %swap the diagonals to put largest in second place
    if vt(3) > vt(2)
        %fprintf('Corner swap, i = %d \n',i)
        vt = vt([1 3 2 4]);
        %vt
    end

    %u(i,:)
    %v(i,:)
    v(i,:) = vt;
end % for

end % function

        




