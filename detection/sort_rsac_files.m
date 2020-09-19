function p=sort_rsac_files(prefix)
% d=sort_rsac_files(file_search)
% example;: log_likelihood('TIFs/','w.mat')
%
% in:
% file_search    directory search string
% d              cell array of file names ordered by time

% insert query to use tifs of Level2 data here
use_tifs = 0;%input_num('Use TIF files? No = 0',0,1);
%obsolete, using only l2 data now
if use_tifs == 1
    d=dir([prefix,'*.tif.mat']);d={d.name};
    if(isempty(d)), error(['No files found for ',prefix]),end
    
    % order the files in time
    nfiles=length(d);
    t=zeros(1,nfiles);
    for i=1:nfiles
        f{i}=[prefix,d{i}];
        t(i)=rsac2time(d{i});
    end
    [t,i]=sort(t);
    p.file={d{i}};
    p.time=t;
    
else
    dhdf=dir([prefix,'*.hdf']);
    dh5 = dir([prefix,'*.h5']);
    dnc = dir([prefix,'*.nc']);
    d=[{dhdf.name},{dh5.name},{dnc.name}];
    if(isempty(d)), error(['No files found for ',prefix]),end
    
    % order the files in time
    nfiles=length(d);
    t=zeros(1,nfiles);
    for i=1:nfiles
        f{i}=[prefix,d{i}];
        t(i)=rsac2time(d{i});
    end
    [t,i]=sort(t);
    p.file={d{i}};
    p.time=t;
    
    %check to make sure geolocation and data files are both present!
    if mod(nfiles,2) ~= 0
        fprintf('Mismatch, number of files is not an even number \n')
    end
    %check that pairs are there delete unpaired granule form list
    match_mask = true(1,nfiles);
    j = 1;
    while j < nfiles-1
        if p.time(j) == p.time(j+1)
            fprintf('time match %s\n',p.file{j})
            fprintf('time match %s\n',p.file{j+1})
            fprintf('Pair %d, %d \n\n',j,j+1)
            j = j+2;
        else
            fprintf('***\n time mismatch %s\n***\n',p.file{j})
            match_mask(j)  = false;
            j = j+1;
        end
    end
    %check for match in last file
    if p.time(nfiles) ~= p.time(nfiles-1)
        fprintf('***\n time mismatch %s\n***\n',p.file{nfiles})
        match_mask(nfiles)  = false;
    end
    %remove unpaired granules
    p.file = p.file(match_mask);
    p.time = p.time(match_mask);
    
    
end % if use_tifs
end %function
