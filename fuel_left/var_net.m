function var_net(net)
%prints out variables needed to make fortran neural network for
%fuel_left_cell_4 subroutine
%clear screen
clc
fprintf('!processing inputs \n\n');
%x1s1off, x1s1gain, x1s1ymin
fort_mat('x1s1off',net.inputs{1}.processSettings{1}.xoffset')
fort_mat('x1s1gain',net.inputs{1}.processSettings{1}.gain')
%fort_mat('x1s1ymin',net.inputs{1}.processSettings{1}.ymin')
fprintf('x1s1ymin = %3.16f \n',net.inputs{1}.processSettings{1}.ymin);

fprintf('\n!layer 1 weights and bias \n')
fort_mat('w1',net.IW{1})
fort_mat('b1',net.b{1})

fprintf('\n!layer 2  weights and bias\n')
fort_mat('w2',net.LW{2,1})
fort_mat('b2',net.b{2})

fprintf('\n!layer 3 weights and bias  \n')
fort_mat('w3',net.LW{3,2})
fort_mat('b3',net.b{3})

fprintf('\n!scale the output \n')
%y1s1ymin, y1s1gain, y1s1off
fprintf('y1s1ymin = %3.16f \n', net.outputs{3}.processSettings{1}.ymin);
fprintf('y1s1gain = %3.16f \n', net.outputs{3}.processSettings{1}.gain);
fprintf('y1s1off =  %3.16f \n', net.outputs{3}.processSettings{1}.xoffset);
fprintf('\n\n\n !done\n')