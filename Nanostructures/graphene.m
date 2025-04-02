% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by select_param to define superlattice parameters specific to graphene (infinite sheet).
% Inputs:
%   r_C = atomic radius of sp2-bonded carbon in Bohr radii
%   a_CC = carbon-carbon bond length in Bohr radii
%   N_z = separation between sheets in sqrt(3)*a_CC (out-of-plane)
%   n_points = number of real- and reciprocal-space points to generate
%   n_ticks = number of independent axis ticks
% Outputs:
%   k = reciprocal-space points at which to compute eigenenergies
%   tau = atomic positions within the supercell
%   r_atom = atomic radii
%   R_gen = generators of superlattice vectors
%   R = lattice translation over which to evaluate wavefunctions
%   plot_struct = container for axis settings
%
% References:
%   [3] M. V. Fischetti and W. G. Vandenberghe. 
%       Advanced Physics of Electron Transport in Semiconductors and Nanostructures.
%       Springer (2016).

function [k, tau, r_atom, R_gen, R, plot_struct] = graphene(r_C, a_CC, N_z, n_points, n_ticks)
    plot_struct = struct();

    % primitive superlattice vectors
    R_gen = sqrt(3)*a_CC*[1/2 -1/2 0; sqrt(3)/2 sqrt(3)/2 0; 0 0 N_z];
    R = (0:1/(n_points-1):1)'*(R_gen(:,1)-R_gen(:,2))';
    plot_struct.R.label = {'$x \ (\sqrt{3}a_{CC})$'};
    plot_struct.R.ticklabels = 0:1/(n_ticks-1):1;
    plot_struct.R.ticks = n_points*plot_struct.R.ticklabels;
    
    % atom positions in primitive cell
    tau = (a_CC/2)*[0 0; -1 1; 0 0];

    % atomic radii
    r_atom = r_C*[1 1];
    
    % irreducible wedge of BZ
    K = 4*pi/(3*sqrt(3)*a_CC);
    M = 2*pi/(3*a_CC);
    D = sqrt(K^2 - M^2);
    k1 = (M:-M/(n_points/3-1):0)'*[sqrt(3)/2 1/2 0];
    nk1 = size(k1,1);
    k2 = (0:K/(n_points/3-1):K)'*[1/2 sqrt(3)/2 0];
    nk2 = size(k2,1);
    k3 = (0:D/(n_points/3-1):D)'*[1/2 -sqrt(3)/2 0];
    nk3 = size(k3,1);
    k = [k1; k2; k3 + ones(size(k3)).*k2(end,:)];
    n_points = size(k,1);
    plot_struct.k.ticks = [1, nk1+1, nk1+nk2+1, n_points-1];
    plot_struct.k.ticklabels = {'$M$' '$\Gamma$' '$K$' '$M$'};
    plot_struct.k.label = {'$k$'};
end