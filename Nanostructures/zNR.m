% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by select_param to define superlattice parameters specific to zigzag nanoribbons (infinite length).
% Inputs:
%   r_C = atomic radius of sp2-bonded carbon in Bohr radii
%   a_CC = carbon-carbon bond length in Bohr radii
%   a_CH = carbon-hydrogen bond length in Bohr radii
%   N_y = separation between ribbons in 3*a_CC (transverse in-plane)
%   N_z = separation between ribbons in 3*a_CC (transverse out-of-plane)
%   N_dw = width of a ribbon in dimer lines
%   n_points = number of real- and reciprocal-space points to generate
%   n_ticks = number of independent axis ticks
% Outputs:
%   k = reciprocal-space points at which to compute eigenenergies
%   tau = atomic positions within the supercell, ordered by carbon first and hydrogen last
%   r_atom = atomic radii, ordered as tau
%   R_gen = generators of superlattice vectors
%   R = lattice translation over which to evaluate wavefunctions
%   plot_struct = container for axis settings
%
% References:
%   [3] M. V. Fischetti and W. G. Vandenberghe. 
%       Advanced Physics of Electron Transport in Semiconductors and Nanostructures.
%       Springer (2016).

function [k, tau, r_atom, R_gen, R, plot_struct] = zNR(r_C, a_CC, a_CH, N_y, N_z, N_dw, n_points, n_ticks)
    plot_struct = struct();

    x_hat = [1; 0; 0]; % unit vector along ribbon length
    y_hat = [0; 1; 0]; % unit vector along ribbon width
    a_hat = [sqrt(3)/2; 1/2; 0]; % unit vector of primitive cell
    N_w = N_dw/2; % number of primitive cells across ribbon width
    n_C = 4*N_w; % number of carbon atoms within supercell
    n_H = 2; % number of hydrogen atoms to terminate dangling edge bonds
    tau = zeros(3,n_C+n_H); 

    % locate carbon atoms within supercell
    tau(:,1) = a_CC*dot(a_hat,x_hat)*x_hat;
    tau(:,2) = a_CC*dot(a_hat,y_hat)*y_hat;
    tau(:,3) = tau(:,2) + a_CC*y_hat;
    tau(:,4) = tau(:,1) + 2*a_CC*y_hat;
    for w = 0:floor(N_w-1) % loop over whole primitive cells along ribbon width
        for n = 1:4 % loop over atoms within primitive cell
            tau(:,w*4+n) = tau(:,n) + a_CC*w*3*y_hat;
        end
    end
    if mod(N_w,1) ~= 0 % odd N_dw; add a half-primitive-cell along ribbon length
        tau(:,4*floor(N_w)+1) = tau(:,1) + floor(N_w)*3*a_CC*y_hat;
        tau(:,4*floor(N_w)+2) = tau(:,2) + floor(N_w)*3*a_CC*y_hat;
    end

    % locate hydrogen atoms along edges within supercell
    tau(:,n_C+1) = tau(:,1) + a_CH*(-y_hat);
    tau(:,n_C+2) = (mod(N_w,1)>0)*(tau(:,2)+floor(N_w)*3*a_CC*y_hat) + (mod(N_w,1)==0)*(tau(:,4)+floor(N_w-1)*3*a_CC*y_hat) + a_CH*y_hat;

    tau = tau - ones(size(tau)).*mean(tau,2); % shift coordinate origin to center of nanostructure

    % atomic radii
    r_atom = ones(1,size(tau,2));
    r_atom(1:n_C) = r_C;

    % generators of superlattice vectors
    R_gen = [sqrt(3)*a_CC 0 0; 0 range(tau(2,:))+N_y*3*a_CC 0; 0 0 N_z*3*a_CC];
    R = (0:1/(n_points-1):1)'*(R_gen(:,1))';
    plot_struct.R.label = {'$x \ (\sqrt{3}a_{CC})$'};
    plot_struct.R.ticklabels = 0:1/(n_ticks-1):1;
    plot_struct.R.ticks = n_points*plot_struct.R.ticklabels;

    % axial k
    k_max = pi/vecnorm(R_gen(:,1));
    dk = k_max/(n_points-1);
    k = (0:dk:k_max)'*[1 0 0];
    plot_struct.k.label = {'$k_x \ (\pi/\sqrt{3}a_{CC})$'};
    plot_struct.k.ticklabels = 0:1/(n_ticks-1):1;
    plot_struct.k.ticks = n_points*plot_struct.k.ticklabels;
end