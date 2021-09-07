function link_namelist(cyc)
%deletes current namelist link and makes new one
if ~exist('cyc','var')
    cyc = input_num('Which namelist to link?',1);
end
ln_str = sprintf('rm namelist.input; ln -s namelist.input_%d namelist.input',cyc);
fprintf('Linking namelist.input_%d\n',cyc);
if system(ln_str)
    printf('Linking failed somehow\n')
end

end