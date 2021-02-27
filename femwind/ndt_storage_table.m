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
        error('symmetric storage not implemented yet')
        for j3=-1:1                  
            for j2=-1:1               
                for j1=-1:1
                    % in fortran won't have the 2+ because
                    % the array t will be indexed -1:1
                    t(1,2+j1,2+j2,2+j3)=0;  % no offset, all here
                    t(2,2+j1,2+j2,2+j3)=0;
                    t(3,2+j1,2+j2,2+j3)=0;
                    % flatten (j1,j2,j3) storage, starting from 1
                    t(4,2+j1,2+j2,2+j3)=1+j1+1+3*(j2+1+3*(j3+1));
                end
            end
        end
    otherwise
        error('unknown storage scheme')
end