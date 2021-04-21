function W=ndt_w_assembly(A,X,u0,lambda,params)

% in: 
%  A  penalty coefficients matrix, size 3x3, s.p.d.
%  X  cell array with 3d arrays x y z coordinates of mesh vertices
%  u0 cell array with 3d arrays initial wind vector at centers in x y z directions
%  lambda scalar field on node; if empty, do not compute K and F
% out:
%  K  stiffness matrix, sparse
%  F  load vector
%  W  gradient of lambda

% do to: now computing all three K,F,W, add a switch to compute only one
% and save time

fprintf('w_assembly:')

n = size(X{1});   % mesh size

% initialize matrices
W = {zeros(n-1),zeros(n-1),zeros(n-1)};

% create matrices
Xloc = zeros(3,8);
kglo=zeros(1,8);
m = n-1;
done_last=0;
disp('Constructing W')
disp('lambda array is')
disp(lambda)

for i3=1:m(3)
    for i2=1:m(2)
        for i1=1:m(1)  % loop over elements
            % now in element (i1,i2,i3)
            for j3=0:1 % loop over corners of the element
                for j2=0:1
                    for j1=0:1   
                        kloc=1+j1+2*(j2+2*j3);  % local index of the node in element (i1, i2. i3)
                        % kloc1=sub2ind([2,2,2],j1+1,j2+1,j3+1);  if kloc1~=kloc, error('kloc'),end 
                        k1 = i1+j1; k2 = i2+j2; k3 = i3+j3; %  position of the node in the global grid
                        kglo(kloc)=k1+n(1)*((k2-1)+n(2)*(k3-1)); % global index
                        % kglo1=sub2ind(n,i1+j1,i2+j2,i3+j3); if kglo1 ~= kglo(kloc), error('kglo'), end
                        for i=1:3
                            Xloc(i,kloc)=X{i}(i1+j1,i2+j2,i3+j3); % node coordinate i
                        end
                    end
                end
            end
            if ~isempty(u0)
                u0loc=[u0{1}(i1,i2,i3),u0{2}(i1,i2,i3),u0{3}(i1,i2,i3)]';
            else
                u0loc=[];
            end
            [~,~,Jg]=hexa(A,Xloc,u0loc); % create local matrix and rhs
            disp('Jg is')
            disp(Jg)
            % instead, accumulate contributions to the global matrix
            if ~isempty(lambda)
                grad = lambda(kglo)'*Jg;  % grad lambda
               %disp('Local lambda is')
               %lambda(kglo)
                %disp(' Grad prior to multiplication by A_inv is')
               %disp(grad) 
                grad = grad/A;
                %disp('Grad after to multiplication by A_inv is')
                 %disp(grad)
                for i=1:3
                    W{i}(i1,i2,i3)=grad(i);
                end
            end
        end
    end
end
if ~isempty(u0)
    for i=1:3
        W{i}=u0{i}+W{i};
    end
end
% fprintf('\n')
end
