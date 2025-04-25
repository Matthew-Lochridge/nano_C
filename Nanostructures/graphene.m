% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by select_param() to define superlattice parameters specific to graphene (infinite sheet).
% Inputs:
%   param = container for nanostructure parameters
%   config = container for figure/axis settings
% Outputs:
%   updated param and config
%
% References:
%   [4] M. V. Fischetti and W. G. Vandenberghe. 
%       Advanced Physics of Electron Transport in Semiconductors and Nanostructures.
%       Springer (2016).

function [param, config] = graphene(param, config)
    
    % atom positions in primitive cell
    param.tau = (param.a_CC/2)*[0 0; -1 1; 0 0];

    % atomic radii
    param.r_atom = param.r_C*[1 1];

    % primitive superlattice vectors
    param.R_gen = sqrt(3)*param.a_CC*[1/2 -1/2 0; sqrt(3)/2 sqrt(3)/2 0; 0 0 param.N_z];

    % real-space path along which to evaluate wavefunctions
    param.R = (0:1/(config.n_points-1):1)'*(param.R_gen(:,1)-param.R_gen(:,2))';
    config.R.label = {'$x \ (\sqrt{3}a_{CC})$'};
    config.R.ticklabels = 0:1/(config.n_ticks-1):1;
    config.R.ticks = config.n_points*config.R.ticklabels;

    K = 4*pi/(3*sqrt(3)*param.a_CC)*[1/2 sqrt(3)/2 0];
    M = 2*pi/(3*param.a_CC)*[0 1 0];
    D = sqrt(vecnorm(K)^2 - vecnorm(M)^2)*[-1 0 0];
    n_step = config.n_points/3 - 1;
    param.dk = vecnorm(D)/n_step;
    
    % reciprocal-space path along boundary of irreducible wedge of BZ
    param.k = [(1:-param.dk/vecnorm(M):1)'*M; (0:param.dk/vecnorm(K):1)'*K; (0:1/n_step:1)'*D + ones(n_step+1,3).*K];
    config.n_points = size(param.k,1);
    config.k.ticks = [1, (n_step+1)+1, 2*(n_step+1)+1, 3*(n_step+1)-1];
    config.k.ticklabels = {'$M$' '$\Gamma$' '$K$' '$M$'};
    config.k.label = {'$k$'};
    %}
    %{
    % full irreducible BZ
    dK = (0:1/n_step:1)'*K;
    IBZ = [];
    for n = 1:n_step+1
        y_slice = (vecnorm(M):-param.dk:dK(n,2))';
        IBZ = [IBZ; [dK(n,1)*ones(size(y_slice)), y_slice, zeros(size(y_slice))]];
    end
    param.k = IBZ;
    config.k.label = {'$k_x$' '$k_y$'};
    config.k.text = {"$\Gamma$" "$M$" "$K$" "$K'$"};
    config.k.sym_points = [zeros(1,3); M; K; K-[vecnorm(K) 0 0]];
    config.vertex = [M; K; K-[vecnorm(K) 0 0]; [vecnorm(K) 0 0]];
    config.sym_proj = [-1 0 0; cos(pi/6) -sin(pi/6) 0; -cos(pi/6) -sin(pi/6) 0; 0 -1 0];
    %}

end