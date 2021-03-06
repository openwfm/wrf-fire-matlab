function DM=mat_wind_flux_div(X,type)
% get matrix of flux divergence by probing its action
% input:
%     X       mesh
%     type    'M': return just the matrix M, 'D': return just the matrix D,
%             'DM': return D*M, default is 'DM'
% output:
%     DM  matrix
if ~exist('type','var')
    type = 'DM';
end

check_mesh(X);
wind_template=grad3z(rand(size(X{1})-1),[1 1 1]);  % cell matrix with entries size of u,v,w
n = sum(cellfun(@numel,wind_template));   % size of vector this will act on

Mfun=@(u)cell2vec(wind2flux(vec2cell(u,wind_template),X)); 
Mmat=fun2mat(Mfun,[n,1,1]);
Dfun=@(u)div3(vec2cell(u,wind_template));
Dmat=fun2mat(Dfun,[n,1,1]);

if (type == 'M')
    DM = Mmat;
elseif (type == 'D')
    DM = Dmat;
else
    DM = Dmat * Mmat;
end

end