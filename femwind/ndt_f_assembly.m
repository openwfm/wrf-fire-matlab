function F=ndt_f_assembly(A,X,u0,params);

% d = size(X,2);    % dimensions
n = size(X{1});   % mesh size in nodes
nn = prod(n);     % total nodes
 

% initialize matrices
F = zeros(nn,1);
Xloc = zeros(3,8);
kglo=zeros(1,8);
uloc = zeros(3);

for ie2=1:n(2)-1  % loop over elements
    for ie3=1:n(3)-1
        for ie1=1:n(1)-1  
            % now in element (ie1,ie2,ie3)
            % build the matrix of node coordinates to pass to hexa
            for ic3=0:1 % loop over corners of the element
                for ic2=0:1
                    for ic1=0:1   
                        iloc=1+ic1+2*(ic2+2*ic3);  % local index of the node in the element
                        k1 = ie1+ic1; k2 = ie2+ic2; k3 = ie3+ic3; %  position of the node in the global grid
                        kglo(iloc)=k1+n(1)*((k2-1)+n(3)*(k3-1)); % global index
                        for i=1:3
                            Xloc(i,iloc)=X{i}(ie1+ic1,ie2+ic2,ie3+ic3); % node coordinate i
                        end
                    end
                end
            end
            if ~isempty(u0)
                u0loc =[u0{1}(ie1,ie2,ie3),u0{2}(ie1,ie2,ie3),u0{3}(ie1,ie2,ie3)]';
            else
                u0loc=[];
            end
            kglo
            [~,Floc,~]=hexa(A,Xloc,u0loc); % compute the local load vector
            for i = 1:8
                F(kglo(i)) = F(kglo(i))+Floc(i);
            end
        end
    end
end
end