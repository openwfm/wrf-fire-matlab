function link_namelist()
%deletes current namelist link and makes new one
cyc = input_num('Which namelist to link?',1);
ln_str = sprintf('rm namelist.input; ln -s namelist.input_%d namelist.input',cyc);
fprintf('Linking namelist.input_%d\n',cyc);
if system(ln_str)
    printf('Linking failed somehow\n')
end

end