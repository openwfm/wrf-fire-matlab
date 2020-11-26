function ts = choose_time_step(f)
%lists time steps in wrfout for user selection
%input   : f - string, path to wrfout file
%output  : ts  - string, timestep

t=nc2struct(f,{'Times'},{});  nframes=size(t.times,2);
alltimes=char(t.times');
for i = 1: nframes
    fprintf('%d :   %s  \n',i,alltimes(i,:))
end
in_t = input_num('Which step to use? ',nframes)
ts = alltimes(in_t,:);

end %function