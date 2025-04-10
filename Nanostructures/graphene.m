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
    
    % reciprocal-space path along boundary of irreducible wedge of BZ
    K = 4*pi/(3*sqrt(3)*param.a_CC);
    M = 2*pi/(3*param.a_CC);
    D = sqrt(K^2 - M^2);
    k1 = (M:-M/(config.n_points/3-1):0)'*[sqrt(3)/2 1/2 0];
    nk1 = size(k1,1);
    k2 = (0:K/(config.n_points/3-1):K)'*[1/2 sqrt(3)/2 0];
    nk2 = size(k2,1);
    k3 = (0:D/(config.n_points/3-1):D)'*[1/2 -sqrt(3)/2 0];
    nk3 = size(k3,1);
    param.k = [k1; k2; k3 + ones(size(k3)).*k2(end,:)];
    config.k.ticks = [1, nk1+1, nk1+nk2+1, nk1+nk2+nk3-1];
    config.k.ticklabels = {'$M$' '$\Gamma$' '$K$' '$M$'};
    config.k.label = {'$k$'};
end