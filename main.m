% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
% 
% Main function to compute energy bands of carbon nanostructures.
% Input:
%   nanostructure = name of nanostructure as a char vector,
%                   e.g. 'Graphene', 'PPP', '5-AGNR', '4-ZGNR', '(5,5)-CNT', etc.
%                   If unspecified, 'graphene' is called by default.
% Output (saved in Data/nanostructure.mat):
%   param = struct containing crystal parameters, energy bands, etc.
%   config = figure/plot settings
%
% References:
%   [1] W. G. Vandenberghe, bulk_pseudo_william.m (UT Dallas)
%   [2] Y. Kurokawa, S. Nomura, T. Takemori, and Y. Aoyagi.
%       "Large-scale calculation of optical dielectric functions of diamond nanocrystallites."
%       Phys. Rev. B 61 12616 (2000).
%   [3] A. Mayer.
%       "Band structure and transport properties of carbon nanotubes using a local pseudopotential and a transfer-matrix technique."
%       Carbon 42 10 (2004).
%   [4] M. V. Fischetti and W. G. Vandenberghe. 
%       Advanced Physics of Electron Transport in Semiconductors and Nanostructures.
%       Springer (2016).
%   [5] J. W. Mintmire and C. T. White. 
%       "Electronic and structural properties of carbon nanotubes."
%       Carbon 33 7 (1995).
%   [6] B. Liu, S. G. Johnson, J. D. Joannopoulos, and L. Lu.
%       "Generalized Gilat–Raubenheimer method for density-of-states calculation in photonic crystals."
%       J. Opt. 20 044005 (2018)

function main(nanostructure)
    addpath('Plot Functions\'); % enable use of functions to plot results
    if nargin == 0
        nanostructure = 'Graphene'; % select graphene by default
    end
    param = struct(); % container for nanostructure parameters
    config = struct(); % container for figure/axis settings

    % constants
    param.Ry = 13.6; % Rydberg energy (eV)
    param.r_H = 5.2917725e-11; % atomic radius of hydrogen (Bohr radius, in m)
    a_d = 3.56683e-10/param.r_H; % diamond lattice constant
    param.r_C = 73e-12/param.r_H; % atomic radius of sp2-bonded carbon
    param.a_CC_graphene = 1.42e-10/param.r_H; % carbon-carbon bond length in graphene
    param.a_CH_graphene = 1.0919e-10/param.r_H; % carbon-hydrogen bond length in methane
    param.a_CC_PPP_hi = 1.478e-10/param.r_H; % carbon-carbon bond length in polyparaphenylene (triply-bonded to carbon)
    param.a_CC_PPP_lo = 1.396e-10/param.r_H; % carbon-carbon bond length in polyparaphenylene (edges)
    param.a_CH_PPP = 1.184e-10/param.r_H; % carbon-hydrogen bond length in polyparaphenylene

    % empirical pseudopotentials normalized to the volume of a carbon atom in diamond [2]
    b_C = [1.781, 1.424, 0.354, 0.938]; % empirical parameters for carbon
    U_C_Kurokawa = @(q_in) (a_d/2)^3*(b_C(1)*(b_C(3)*q_in.^2-b_C(2))./(exp(b_C(3)*q_in.^2-b_C(4))+1)); % pseudopotential function for carbon
    b_H = [-0.397, 0.0275, 0.1745, -0.0531, 0.0811, -1.086, 2.71, -2.86]; % empirical parameters for hydrogen
    U_H_Kurokawa = @(q_in) (a_d/2)^3*((q_in<=2).*(b_H(1) + b_H(2)*q_in + b_H(3)*q_in.^2 + b_H(4)*q_in.^3) + (q_in>2).*(b_H(5)./(q_in+(q_in==0)) + b_H(6)./(q_in.^2+(q_in==0)) + b_H(7)./(q_in.^3+(q_in==0)) + b_H(8)./(q_in.^4+(q_in==0)))); % pseudopotential function for hydrogen
    
    % empirical pseudopotentials [3] normalized to the volume of a carbon atom in diamond [4]
    A_j = [0.7796, 2.1837, -7.2698]; % empirical parameters for carbon
    a_j = [0.12126, 1.9148, 0.60078]; % empirical parameters for carbon
    b_j = [-2.68574, -0.11989, 2.27101]; % -A_j.*(pi./a_j).^(3/2); % Fourier coefficients
    U_C_Mayer = @(q_in) (a_d/2)^3*(b_j(1)*exp(-q_in.^2/(4*a_j(1))) + b_j(2)*exp(-q_in.^2/(4*a_j(2))) + b_j(3)*exp(-q_in.^2/(4*a_j(3)))); % pseudopotential function for carbon

    % cutoff parameters
    param.E_cut = 15; % cutoff energy (Ry)
    param.max_G = 100; % maximum manhattan distance of G vectors

    % supercell size parameters
    param.N_a = 1; % number of primitive motifs along a tube within the supercell (axial); Set to 1 for an infinite tube.
    param.N_x = 0; % separation between finite ribbons or tubes in primitive translations (axial); Set to 0 for an infinite ribbon or tube.
    param.N_y = 4; % separation between ribbons in primitive translations (transverse in-plane)
    param.N_z = 15e-10/param.r_H; % separation between graphene sheets and ribbons in primitive translations (transverse out-of-plane), or supercell cross-section length for tubes

    % figure/axis settings
    param.nanostructure = nanostructure;
    config.E.lim = 5; % maximum energy shown in band plot (eV)
    config.n_points = 101; % number of real- and reciprocal-space points
    config.n_ticks = 6;
    config.interpreter = 'latex';
    set(groot,'defaultTextInterpreter',config.interpreter);
    set(groot,'defaultAxesTickLabelInterpreter',config.interpreter);

    % determine parameters corresponding to specific nanostructures
    [param, config] = select_param(param, config);
    if ~isfield(param,'tau')
        disp('Quasicrystal parameters cannot be determined.')
        return;
    end
    if ~strcmpi(nanostructure,'graphene') % display supercell
        plot_supercell(param);
    end

    % compute energy bands and corresponding wavefunctions
    param = bands(param, U_C_Kurokawa, U_H_Kurokawa);
    
    % save data
    disp('Saving data...');
    tic
    save(append('Data/',nanostructure),'-v7.3');
    toc

    % plot results
    disp('Plotting band structure...');
    tic
    plot_bands(param, config);
    toc
    disp('Plotting wavefunctions...')
    tic
    plot_wavefunc(param);
    toc
    disp('Finished.');
end