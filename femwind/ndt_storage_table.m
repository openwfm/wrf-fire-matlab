function t=ndt_storage_table(m)
% t=ndt_storage_table(m)
% returns 6 storage indices for each of 27 neigbors as array size (6,3,3,3)
% t(1:3,...) row offsets 0 = this this row 

switch m
    case 27
        % no row offset
        t=zeros(4,3,3,3);
        % flat indexing within the row
        t(4,:,:,:)=reshape(1:27,3,3,3);
    case 14
        g=@(j1,j2,j3)1+(j1+1)+3*(j2+1+3*(j3+1)); % storage function
        g0 = g(0,0,0);
        n=0;
        for j3=-1:1                  
            for j2=-1:1               
                for j1=-1:1
                    gj=g(j1,j2,j3);
                    if gj >= g0
                        % we are in the upper triangle including diagonal
                        % connecting from 0 0 0 to j1 j2 j3
                        n=n+1;
                        t(1:3,2+j1,2+j2,2+j3)=0; % is centered at 2 not 0
                        t(4  ,2+j1,2+j2,2+j3)=n;
                        % disp([gj,j1,j2,j3,squeeze(t(:,2+j1,2+j2,2+j3))'])
                    end
                end
            end
        end
        if n ~= 14  % just checking
            nu, error('wrong size of upper triangular + diagonal part')
        end
        % now fill the lower triangular entries by index  
        % to the upper triangular entry on the other row
        for j3=-1:1                  
            for j2=-1:1               
                for j1=-1:1
                    gj=g(j1,j2,j3);
                    if gj < g0
                        % we are in lower triangle
                        % conecting from j1 j2 j3 to 0 0 0
                        t(1:3,2+j1,2+j2,2+j3)=[j1; j2; j3]; % from
                        t(4  ,2+j1,2+j2,2+j3)=t(4  ,2-j1,2-j2,2-j3); % to
                        % disp([gj,j1,j2,j3,squeeze(t(:,2+j1,2+j2,2+j3))'])
                        if t(4  ,2+j1,2+j2,2+j3) == 0
                            error('no target')
                        end
                    end
                end
            end
        end
    otherwise
        m,error('unknown storage scheme, must be 14 or 27')
end