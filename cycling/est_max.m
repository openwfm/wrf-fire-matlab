function perc = est_max(ps,r)
%estimate best cut-off for ROS when fitting model
%mdl2_1 = fitlm(pc17.red.fmc_g(msk_1),r2(msk_1))
%mdl.Coefficients(1,4),(2,4) are the pvaleus

stopping = 0;
for i = 100:-1:1;%1:100
r_cut(i) = i/500;
msk1 = r<r_cut(i);
msk2 = r > 0.0001;
msk = logical(msk1.*msk2);
s(i) = sum(msk(:));
if s(i) > 20
    mdl = fitlm(ps.red.fmc_g(msk),r(msk));
    p(i)=table2array(mdl.Coefficients(2,4));
    if p(i) <= 0.01 && stopping == 0
        fprintf('cutoff at ros = %f \n',r_cut(i))
        stopping = 1;
        perc = r_cut(i);
        break
    end
end

end%for i
% figure,plot(r_cut,p)
% title('P-values')
% ylim([0 0.1])
% 
% 
% figure,plot(r_cut,s)
% title('samples')
if perc == []
    perc = 0.1
end

end%function