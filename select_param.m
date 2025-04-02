% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
% 
% Function called by nano_C to select which structural parameters to generate.
% Input:
%   r_C = atomic radius of sp2-bonded carbon in Bohr radii
%   a_CC = carbon-carbon bond length in Bohr radii
%   a_CH = carbon-hydrogen bond length in Bohr radii
%   N_a = number of primitive cells along a tube within the supercell (axial); Input 1 for an infinite tube.
%   N_x = separation between finite ribbons or tubes in a_CC (axial); Input 0 for an infinite tube.
%   N_y = separation between ribbons or tubes in a_CC (transverse, in-plane for ribbons)
%   N_z = separation between ribbons, tubes, or infinite graphene sheets in a_CC (transverse, out-of-plane for ribbons and graphene)
%   structure = name of selected nanostructure as a char vector,
%               e.g. '9-AGNR', '8-ZGNR', '(10,10)@(15,15)-CNT', etc.
%               If unspecified, 'graphene' is called by default.
%   n_points = number of real- and reciprocal-space points
%   n_ticks = number of independent axis ticks
% Output:
%   k = reciprocal-space points at which eigenenergies are computed
%   tau = atomic positions within the supercell
%   r_atom = atomic radii ordered as tau
%   R_gen = primitive superlattice vectors
%   R = lattice translation over which to evaluate wavefunctions
%   plot_struct = container for axis settings

function [k, tau, r_atom, R_gen, R, plot_struct] = select_param(r_C, a_CC, a_CH, N_a, N_x, N_y, N_z, structure, n_points, n_ticks)
    addpath('Nanostructures\'); % enable use of functions to generate nanostructures

    [allotrope, N_dl, N_dw, N_1, N_2] = parse_structure(structure); % decompose structure into type and parameters

    % define k, tau, r_atom, R_gen, and plot settings based on type and parameters
    if strcmpi(allotrope,'graphene')
        disp('Graphene selected (default).')
        [k, tau, r_atom, R_gen, R, plot_struct] = graphene(r_C, a_CC, N_z, n_points, n_ticks);

    elseif (N_dw > 0) && isscalar(N_dl) % nanoribbons
        disp(append(structure,' selected.'))
        if (N_dl == 2) && strcmpi(allotrope,'agnr')
            [k, tau, r_atom, R_gen, R, plot_struct] = aNR(r_C, a_CC, a_CH, N_y, N_z, N_dw, n_points, n_ticks);
        elseif (N_dl == 2) && strcmpi(allotrope,'zgnr')
            [k, tau, r_atom, R_gen, R, plot_struct] = zNR(r_C, a_CC, a_CH, N_y, N_z, N_dw, n_points, n_ticks);
        elseif strcmpi(allotrope,'gnr')
            [k, tau, r_atom, R_gen, R, plot_struct] = nanoribbon(r_C, a_CC, a_CH, N_x, N_y, N_z, N_dl, N_dw, n_points, n_ticks);
        end

    elseif length(N_1)==length(N_2) && sum(N_1+N_2>0)==length(N_1) % nanotubes
        disp(append(structure,' selected.'))
        k = [];
        tau = [];
        r_atom = [];
        R_gen = [];
        R = [];
        plot_struct = [];
        circ = zeros(size(N_1));
        circ_func = @(n,m) sqrt(n^2+m^2+n*m); % circumference of (n,m)-CNT in sqrt(3)*a_CC
        for t = 1:length(N_1) % loop over coaxial tubes
            if N_1(t) == 0 || N_2(t) == 0 % zigzag
                [k_t, tau_t, r_atom_t, R_gen_t, R_t, plot_struct_t] = zNT(r_C, a_CC, a_CH, N_a, N_x, N_z, max(abs([N_1(t) N_2(t)])), n_points, n_ticks);
            elseif N_2(t) == N_1(t) % armchair
                [k_t, tau_t, r_atom_t, R_gen_t, R_t, plot_struct_t] = aNT(r_C, a_CC, a_CH, N_a, N_x, N_z, N_1(t), n_points, n_ticks);
            else % chiral
                [k_t, tau_t, r_atom_t, R_gen_t, R_t, plot_struct_t] = nanotube(r_C, a_CC, a_CH, N_a, N_x, N_z, N_1(t), N_2(t), n_points, n_ticks);
            end
            % append atoms to supercell for each tube
            tau = [tau tau_t];
            r_atom = [r_atom r_atom_t];
            % use only k, R_gen, and R for the largest tube
            circ(t) = circ_func(N_1(t),N_2(t));
            if circ(t) == max(circ)
                k = k_t;
                R_gen = R_gen_t;
                R = R_t;
                plot_struct = plot_struct_t;
            end
        end

    else
        k = [];
        tau = [];
        r_atom = [];
        R_gen = [];
        R = [];
        plot_struct = [];
        return;
    end
end

