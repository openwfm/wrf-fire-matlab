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
P = p.p+p.pb;            % air pressure
theta = p.t+300;             % potential temperature in F
% theta = T*(p0/p)^c, p0=1e5, c=0.286 https://en.wikipedia.org/wiki/Potential_temperature    
temp = theta .* (P*1e-5).^0.286;
TV=(1+0.61*p.qvapor).*temp; % virtual temperature
RHO=P./(287*TV);           % air density
EL =(p.ph + p.phb)/9.81; % elevation at w-points
Z = (EL(:,:,1:end-1)+EL(:,:,2:end))/2; % elevation at centers

dp=mean(P(:,:,1)-p.psfc,'all');
sp=mean(p.psfc,'all');
fprintf('average difference P(:,:,1)-PSFC = %g = %g %s\n',dp,100*dp/sp,'%')
dt=mean(temp(:,:,1)-p.t2,'all');
mt2=big(p.t2);
fprintf('average difference temp(:,:1) - T2 = %g max abs t2=%g\n',dt,mt2)
end
