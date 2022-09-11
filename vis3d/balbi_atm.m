function q=balbi_atm(f,s)
% q=balbi_atm(f,s)
% 
% Get atmospheric inputs for the Balbi model from a wrfout
% 
% in:
%   f   wrfout file namee
%   s   step (frame number in the file)
% out:
%   q   structure
% 
% based on formulas from Derek

if nargin<2
    s=1;
end
p=nc2struct(f,{'QVAPOR','T','P','PB','T2','PSFC','PH','PHB'},{},s);
T = p.t+300;             % temperature in F
TV=(1+0.61*p.qvapor).*T; % virtual temperature
P = p.p+p.pb;            % air pressure
RHO=P./(287*TV);           % air density
EL =(p.ph + p.phb)/9.81; % elevation at w-points
Z = (EL(:,:,1:end-1)+EL(:,:,2:end))/2; % elevation at centers

dp=mean(P(:,:,1)-p.psfc,'all');
sp=mean(p.psfc,'all');
fprintf('average difference P(:,:,1)-PSFC = %g = %g %s\n',dp,100*dp/sp,'%')
dt=mean(T(:,:,1)-p.t2,'all');
mt2=big(p.t2);
fprintf('average difference T(:,:1)+300 - T2 = %g max abs t2=%g\n',dt,mt2)
end
