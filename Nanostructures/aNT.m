% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by select_param to define superlattice parameters specific to armchair nanotubes.
% Inputs:
%   r_C = atomic radius of sp2-bonded carbon in Bohr radii
%   a_CC = carbon-carbon bond length in Bohr radii
%   a_CH = carbon-hydrogen bond length in Bohr radii
%   N_a = number of primitive motifs along a tube within the supercell (axial); Input 1 for an infinite tube.
%   N_x = separation between finite tubes in sqrt(3)*a_CC (axial); Input 0 for an infinite tube.
%   N_z = separation between tubes in 3*a_CC (transverse)
%   N_1 = number of each primitive graphene translations making up the wrapping vector
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

function [k, tau, r_atom, R_gen, R, plot_struct] = aNT(r_C, a_CC, a_CH, N_a, N_x, N_z, N_1, n_points, n_ticks)
    plot_struct = struct();

    x_hat = [1; 0; 0]; % unit vector along tube length (axial)
    r_hat = @(phi) [0; cos(phi); sin(phi)]; % radial unit vector as a function of azimuthal angle phi
    R_graphene = sqrt(3)*a_CC*[1/2 -1/2 0; sqrt(3)/2 sqrt(3)/2 0]; % primitive lattice vectors of graphene
    C = N_1*(R_graphene(:,1) + R_graphene(:,2)); % wrapping vector (circumference)
    r_tube = norm(C)/(2*pi); % nanotube radius
    theta = a_CC/r_tube; % characteristic angle of carbon-carbon bond
    n_C = 4*N_1*N_a; % number of carbon atoms within supercell
    n_H = 4*N_1*(N_x>0); % number of hydrogen atoms to terminate dangling bonds
    tau = zeros(3,n_C+n_H);
    
    % locate carbon atoms within supercell
    idx_C = 0; % number of carbon atoms counted so far
    for l = 0:floor(N_a-1) % loop over tube length within supercell
        for c = 0:N_1-1 % loop over tube circumference
            tau(:,idx_C+1) = r_tube*r_hat(c*3*theta) + l*sqrt(3)*a_CC*x_hat;
            tau(:,idx_C+2) = r_tube*r_hat((c*3+1/2)*theta) + (l*sqrt(3)+sqrt(3)/2)*a_CC*x_hat;
            tau(:,idx_C+3) = r_tube*r_hat((c*3+3/2)*theta) + (l*sqrt(3)+sqrt(3)/2)*a_CC*x_hat;
            tau(:,idx_C+4) = r_tube*r_hat((c*3+2)*theta) + l*sqrt(3)*a_CC*x_hat;
            idx_C = idx_C + 4;
        end
    end
    if mod(N_a,1) > 0 % half-integer N_a; add half-motif along circumference
        for c = 0:N_1-1
            tau(:,idx_C+1) = r_tube*r_hat(c*3*theta) + floor(N_a)*sqrt(3)*a_CC*x_hat;
            tau(:,idx_C+2) = r_tube*r_hat((c*3+2)*theta) + floor(N_a)*sqrt(3)*a_CC*x_hat;
            idx_C = idx_C + 2;
        end
    end

    if n_H > 0 % finite tubes; locate hydrogen atoms along edges within supercell
        idx_H = 0; % number of hydrogen atoms counted so far
        for c = 0:N_1-1
            tau(:,n_C+idx_H+1) = tau(:,4*c+1) + a_CH*(tau(:,4*c+2)-tau(:,4*c+1))/norm(tau(:,4*c+2)-tau(:,4*c+1)).*[-1; 1; 1];
            tau(:,n_C+idx_H+2) = tau(:,4*c+4) + a_CH*(tau(:,4*c+3)-tau(:,4*c+4))/norm(tau(:,4*c+3)-tau(:,4*c+4)).*[-1; 1; 1];
            tau(:,n_C+idx_H+3) = (mod(N_a,1)>0)*(tau(:,4*c+1)+floor(N_a)*sqrt(3)*a_CC*x_hat+a_CH*(tau(:,4*c+2)-tau(:,4*c+1))/norm(tau(:,4*c+2)-tau(:,4*c+1))) + (mod(N_a,1)==0)*(tau(:,4*c+2)+floor(N_a-1)*sqrt(3)*a_CC*x_hat+a_CH*(tau(:,4*c+1)-tau(:,4*c+2))/norm(tau(:,4*c+1)-tau(:,4*c+2)).*[-1; 1; 1]);
            tau(:,n_C+idx_H+4) = (mod(N_a,1)>0)*(tau(:,4*c+4)+floor(N_a)*sqrt(3)*a_CC*x_hat+a_CH*(tau(:,4*c+3)-tau(:,4*c+4))/norm(tau(:,4*c+3)-tau(:,4*c+4))) + (mod(N_a,1)==0)*(tau(:,4*c+3)+floor(N_a-1)*sqrt(3)*a_CC*x_hat+a_CH*(tau(:,4*c+4)-tau(:,4*c+3))/norm(tau(:,4*c+4)-tau(:,4*c+3)).*[-1; 1; 1]);
            idx_H = idx_H + 4;
        end
    end

    tau = tau - ones(size(tau)).*mean(tau,2); % shift coordinate origin to center of nanostructure

    % atomic radii
    r_atom = ones(1,size(tau,2));
    r_atom(1:n_C) = r_C;

    % generators of superlattice vectors
    tau_range = range(tau,2);
    R_gen = [(N_x>0)*(tau_range(1)+N_x*sqrt(3)*a_CC)+(N_x==0)*sqrt(3)*a_CC 0 0; 0 tau_range(2)+N_z*3*a_CC 0; 0 0 tau_range(3)+N_z*3*a_CC];
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