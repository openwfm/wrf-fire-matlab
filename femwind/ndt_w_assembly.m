function F=ndt_w_assembly(A,X,u0,lambda,params);

% d = size(X,2);    % dimensions
n = size(X{1});   % mesh size in nodes
nn = prod(n);     % total nodes
 

% initialize matrices
W = {zeros(n-1),zeros(n-1),zeros(n-1)};
Xloc = zeros(3,8);
kglo=zeros(1,8);

for ie3=1:n(3)-1  % loop over elements
    for ie2=1:n(2)-1
        for ie1=1:n(1)-1  
            % now in element (ie1,ie2,ie3)
            % build the matrix of node coordinates to pass to hexa
            for ic3=0:1 % loop over corners of the element
                for ic2=0:1
                    for ic1=0:1   
                        iloc=1+ic1+2*(ic2+2*ic3);  % local index of the node in the element
                        k1 = ie1+ic1; k2 = ie2+ic2; k3 = ie3+ic3; %  position of the node in the global grid
                        kglo(iloc)=k1+n(1)*((k2-1)+n(2)*(k3-1)); % global index
                        for i=1:3
                            Xloc(i,iloc)=X{i}(ie1+ic1,ie2+ic2,ie3+ic3); % node coordinate i
                        end
                    end
                end
            end
            [~,~,Jg]=hexa(A,Xloc,u0); % compute the local load vector
            
            grad = ;
            grad = ;
            for i=1:3
                    W{i}(i1,i2,i3)=grad(i);
            end
        end
    end
end
end
