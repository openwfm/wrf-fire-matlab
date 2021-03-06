function [U,varargout]=mass_cons_flux(U0,X,method,check)
% mass consistent flux approximation
% input:
%   U0         cell array, wind vectors u v w on staggered grids
%   X          cell array, nodal coordinatees
%   method     'direct' or 'pcg' method for solving system
%   check      flag for optional timing and error checking
% output:
%   U          cell array, wind vectors u v w
%   Lambda     pressure field (optional)
%   
%
% using:
%   D = divergence of flux on cell sides
%   M = transformation of wind to flux through cell sides with
%       zero flux boundary condition at the bottom (the terrain)
%   A = weights of the wind vector components, U'*A*U approx ||U||_L2^2
% 
% compute:
%
%    min 1/2 (U - U0)'*A*(U - U0) s.t. D*M*U = 0
%
%    <=>    A*(U - U0) + M'*D'*Lambda = 0
%
%           D*M*U                     = 0
% 
%    <=>    U = U0 - inv(A)*M'*D'*Lambda
%
%           D*M*inv(A)*M'*D'*Lambda = D*M*U0

if ~exist('check','var')
    check = false;
end

if check
    % Timing start
    tstart=tic;
end

% Creates the weight vector... diagonal for now?
hvr = 1e4; % ratio of horizontal and vertical penalty: Sherman 1978 page 315
s = cell_sizes(X);
Avec = cell2vec({hvr*s.weight_u,hvr*s.weight_v,s.weight_w});
Ainv = diag(sparse(1./Avec));

%method = 'iterative'
switch method
    case 'direct'
        DM = mat_gen_wind_flux_div(X); % divergence of u0
        U0vec = cell2vec(U0);
        rhs = DM * U0vec;
        L = DM * Ainv * DM';
        Lambda = L \ rhs;
        Uvec = U0vec - Ainv * DM' * Lambda;
        U = vec2cell(Uvec,U0);
    case 'pcg'
        % uses pcg and implicit multiplication
        DM = mat_gen_wind_flux_div(X); % divergence of u0
        U0vec = cell2vec(U0);
%         rhs = Mmul_v(X, 'n', U0vec);
%         rhs = Dmul_v(X, 'n', rhs);
        rhs = DM * U0vec;
        L = DM * Ainv * DM';
%        PC_L = ichol(PC,struct('michol','on'));
%        lhs_mat = @(x) lhs_apply(X, Ainv, x);
%        Lambda = pcg(lhs_mat, rhs, 1e-10, 100, PC_L, PC_L');
%        Lambda = pcg(lhs_mat, rhs, 1e-8, 400);
%        [Lambda,FLAG,RELRES,ITER,RESVEC] = pcg(lhs_mat, rhs, 1e-6, 400);
        [Lambda,FLAG,RELRES,ITER,RESVEC] = pcg(L, rhs, 1e-6, 400);
        figure;
        semilogy(RESVEC)
        title('Residuals using PCG')
        xlabel('Iterations'),ylabel('Residuals')
        Uvec = U0vec - Ainv * DM' * Lambda;
        U = vec2cell(Uvec,U0);
    otherwise
        error('only direct and iterative methods implemented')
end 

if check
    % Check that the divergence of the computed U is 0
    err_div = big(div3(wind2flux(U,X)))
    varargout{2} = err_div;
    fprintf('mass_cons_flux %g seconds\n',toc(tstart))
end
if length(varargout) > 0
    varargout{1} = Lambda;
end
% check no correction in normal speed at the bottom
% if check, err_corr_w_bottom = big(u{3}(:,:,1)-U0{3}(:,:,1)), end
end
