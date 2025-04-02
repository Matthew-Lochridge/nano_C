% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
% 
% Main function to compute energy bands of carbon nanostructures.
% Input:
%   structure = name of nanostructure as a char vector,
%               e.g. '9-AGNR', '8-ZGNR', '(5,5)-CNT', etc.
%               If unspecified, 'graphene' is called by default.
% Output:
%   E = array of energy bands in Ry ordered as (k,G)
%   psi = cell of wavefunctions ordered by increasing eigenenergy
%
% References:
%   [1] W. G. Vandenberghe. 
%       bulk_pseudo_william.m.
%   [2] Y. Kurokawa, S. Nomura, T. Takemori, and Y. Aoyagi.
%       "Large-scale calculation of optical dielectric functions of diamond nanocrystallites."
%       Phys. Rev. B 61 12616 (2000).
%   [3] M. V. Fischetti and W. G. Vandenberghe. 
%       Advanced Physics of Electron Transport in Semiconductors and Nanostructures.
%       Springer (2016).
%   [4] J. W. Mintmire and C. T. White. 
%       "Electronic and structural properties of carbon nanotubes."
%       Carbon 33 7 (1995).
%   [5] J. Fang, W. G. Vandenberghe, and M. V. Fischetti. 
%       "Microscopic dielectric permittivities of graphene nanoribbons and graphene."
%       Phys. Rev. B 94 045318 (2016).

function [E, psi] = nano_C(structure)
    addpath('Plot Functions\'); % enable use of functions to plot results
    if nargin == 0
        structure = 'Graphene'; % select graphene by default
    end

    % constants
    Ry = 13.6; % Rydberg energy (eV)
    r_H = 5.2917725e-11; % atomic radius of hydrogen (Bohr radius, in m)
    r_C = 73e-12/r_H; % atomic radius of sp2-bonded carbon
    a_CC = 1.42e-10/r_H; % carbon-carbon bond length in graphene
    a_CH = 1.0919e-10/r_H; % carbon-hydrogen bond length in methane

    % empirical pseudopotential parameters
    E_cut = 15; % cutoff energy (Ry)
    max_G = 20; % maximum manhattan distance of G vectors

    % supercell size parameters
    N_a = 1; % number of primitive motifs along a tube within the supercell (axial); Set to 1 for an infinite tube.
    N_x = 0; % separation between finite ribbons or tubes in primitive translations (axial); Set to 0 for an infinite ribbon or tube.
    N_y = 4; % separation between ribbons in primitive translations (transverse in-plane)
    N_z = 5; % separation between ribbons, tubes, or graphene sheets in primitive translations (transverse, out-of-plane for ribbons and graphene)

    % plot settings
    E_range = 20; % maximum energy shown in band plot (eV)
    n_points = 300; % number of real- and reciprocal-space points
    n_ticks = 6; % number of independent axis ticks

    % determine parameters corresponding to structure
    [k, tau, r_atom, R_gen, R, plot_struct] = select_param(r_C, a_CC, a_CH, N_a, N_x, N_y, N_z, structure, n_points, n_ticks);
    if isempty(tau)
        disp('Quasicrystal parameters cannot be determined.')
        E = [];
        psi = [];
        return;
    end
    if ~strcmpi(structure,'graphene') % display supercell
        plot_supercell(structure, 10^9*r_H*a_CC, 10^9*r_H*a_CH, 10^9*r_H*tau, r_atom, 10^9*r_H*a_CC*N_x, 10^9*r_H*a_CC*N_y, 10^9*r_H*a_CC*N_z);
    end

    % compute energy bands and corresponding wavefunctions
    [E, psi] = bands(real(k), R_gen, R, tau, r_atom, E_cut, max_G);
    plot_bands(E*Ry, E_range, structure, plot_struct);
    plot_wavefunc(psi, structure, plot_struct);
end