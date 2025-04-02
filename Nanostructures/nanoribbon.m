% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by select_param to define superlattice parameters specific to finite nanoribbons.
% Inputs:
%   r_C = atomic radius of sp2-bonded carbon in Bohr radii
%   a_CC = carbon-carbon bond length in Bohr radii
%   a_CH = carbon-hydrogen bond length in Bohr radii
%   N_x = separation between finite ribbons in 3*a_CC (axial)
%   N_y = separation between ribbons in sqrt(3)*a_CC (transverse in-plane)
%   N_z = separation between ribbons in sqrt(3)*a_CC (transverse out-of-plane)
%   N_dl, N_dw = length (axial) and width (transverse in-plane) of a ribbon within the supercell in dimer lines
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

function [k, tau, r_atom, R_gen, R, plot_struct] = nanoribbon(r_C, a_CC, a_CH, N_x, N_y, N_z, N_dl, N_dw, n_points, n_ticks)
    plot_struct = struct();
    
    x_hat = [1; 0; 0]; % unit vector along ribbon length
    y_hat = [0; 1; 0]; % unit vector along ribbon width
    a_hat = [1/2; sqrt(3)/2; 0]; % unit vector of primitive cell
    N_l = N_dl/2; % number of primitive cells across ribbon length within the supercell
    N_w = N_dw/2; % number of primitive cells across ribbon width within the supercell
    n_C = 4*N_l*N_w; % number of carbon atoms within supercell
    n_H = (N_x>0)*2*N_w + (N_y>0)*4*N_l; % number of hydrogen atoms to terminate dangling edge bonds
    tau = zeros(3,n_C+n_H);
    
    % locate carbon atoms within supercell
    idx_C = 0; % number of carbon atoms counted so far
    tau(:,2) = a_CC*a_hat;
    tau(:,3) = tau(:,2) + a_CC*x_hat;
    tau(:,4) = tau(:,1) + 2*a_CC*x_hat;
    for l = 0:floor(N_l-1) % loop over whole primitive cells along ribbon length, if finite
        for w = 0:floor(N_w-1) % loop over whole primitive cells along ribbon width, if finite
            for n = 1:4 % loop over atoms within primitive cell
                tau(:,idx_C+n) = tau(:,n) + a_CC*(l*3*x_hat + w*sqrt(3)*y_hat);
            end
            idx_C = idx_C + 4;
        end
    end
    if mod(N_l,1) > 0 % odd N_dl; add floor(N_w) half-primitive-cells along ribbon width
        for w = 0:floor(N_w-1)
            tau(:,idx_C+1) = tau(:,1) + a_CC*(floor(N_l)*3*x_hat + w*sqrt(3)*y_hat);
            tau(:,idx_C+2) = tau(:,2) + a_CC*(floor(N_l)*3*x_hat + w*sqrt(3)*y_hat);
            idx_C = idx_C + 2;
        end
    end
    if mod(N_w,1) > 0 % odd N_dw; add floor(N_l) half-primitive-cells along ribbon length
        for l = 0:floor(N_l-1)
            tau(:,idx_C+1) = tau(:,1) + a_CC*(l*3*x_hat + floor(N_w)*sqrt(3)*y_hat);
            tau(:,idx_C+2) = tau(:,4) + a_CC*(l*3*x_hat + floor(N_w)*sqrt(3)*y_hat);
            idx_C = idx_C + 2;
        end
        if mod(N_l,1) > 0 % odd N_dl; add one more atom to ribbon corner
            tau(:,idx_C+1) = tau(:,1) + a_CC*(floor(N_l)*3*x_hat + floor(N_w)*sqrt(3)*y_hat);
            idx_C = idx_C + 1;
        end
    end

    if n_H > 0 % locate hydrogen atoms along ribbon edges
        idx_H = 0; % number of hydrogen atoms counted so far
        if N_x > 0 % finite length
            for w = 0:floor(N_w-1)
                tau(:,n_C+idx_H+1) = tau(:,1) + w*sqrt(3)*a_CC*y_hat - a_CH*x_hat;
                tau(:,n_C+idx_H+2) = (mod(N_l,1)>0)*(tau(:,2)+a_CC*(floor(N_l)*3*x_hat+w*sqrt(3)*y_hat)) + (mod(N_l,1)==0)*(tau(:,4)+a_CC*(floor(N_l-1)*3*x_hat+w*sqrt(3)*y_hat)) + a_CH*x_hat;
                idx_H = idx_H + 2;
            end
            if mod(N_w,1) > 0 % odd N_dw; add one atom to a ribbon corner
                tau(:,n_C+idx_H+1) = tau(:,1) + a_CC*floor(N_w)*sqrt(3)*y_hat + a_CH*(-x_hat);
                idx_H = idx_H + 1;
            end
        end
        if N_y > 0 % finite width
            for l = 0:floor(N_l-1)
                tau(:,n_C+idx_H+1) = tau(:,1) + l*3*a_CC*x_hat + a_CH*(x_hat-a_hat);
                tau(:,n_C+idx_H+2) = tau(:,4) + l*3*a_CC*x_hat + a_CH*(-a_hat);
                tau(:,n_C+idx_H+3) = (mod(N_w,1)>0)*(tau(:,1)+a_CC*(l*3*x_hat+floor(N_w)*sqrt(3)*y_hat)+a_CH*a_hat) + (mod(N_w,1)==0)*(tau(:,2)+a_CC*(l*3*x_hat+floor(N_w-1)*sqrt(3)*y_hat)+a_CH*(a_hat-x_hat));
                tau(:,n_C+idx_H+4) = (mod(N_w,1)>0)*(tau(:,4)+a_CC*(l*3*x_hat+floor(N_w)*sqrt(3)*y_hat)+a_CH*(a_hat-x_hat)) + (mod(N_w,1)==0)*(tau(:,3)+a_CC*(l*3*x_hat+floor(N_w-1)*sqrt(3)*y_hat)+a_CH*a_hat);
                idx_H = idx_H + 4;
            end
            if mod(N_l,1) > 0 % odd N_dl; add one atom to a ribbon corner
                tau(:,n_C+idx_H+1) = tau(:,1) + a_CC*floor(N_l)*3*x_hat + a_CH*(x_hat-a_hat);
                idx_H = idx_H + 1;
            end
        end
        if (N_x>0) && (N_y>0) && (mod(N_l,1)>0) && (mod(N_w,1)>0) % add one more atom to a ribbon corner
            tau(:,n_C+idx_H+1) = tau(:,1) + a_CC*(floor(N_l)*3*x_hat + floor(N_w)*sqrt(3)*y_hat) + a_CH*a_hat;
            idx_H = idx_H + 1;
        end
    end

    tau = tau - ones(size(tau)).*mean(tau,2); % shift coordinate origin to center of nanostructure

    % atomic radii
    r_atom = ones(1,size(tau,2));
    r_atom(1:n_C) = r_C;

    % generators of superlattice vectors
    tau_range = range(tau,2);
    R_gen = [tau_range(1)+N_x*3*a_CC 0 0; 0 tau_range(2)+N_y*sqrt(3)*a_CC 0; 0 0 N_z*sqrt(3)*a_CC];
    R = (0:1/(n_points-1):1)'*(R_gen(:,1) + 1i*R_gen(:,2))';
    plot_struct.R.label = {'$x \ (3a_{CC})$','$y \ (\sqrt{3}a_{CC})$'};
    plot_struct.R.ticklabels = 0:1/(n_ticks-1):1;
    plot_struct.R.ticks = n_points*plot_struct.R.ticklabels;
        
    % k_x (armchair, Re) and k_y (zigzag, Im)
    k_max = pi./[vecnorm(R_gen(:,1)); vecnorm(R_gen(:,2))];
    dk = k_max/(n_points-1);
    k = (0:dk(1):k_max(1))'*[1 0 0] + (0:dk(2):k_max(2))'*[0 1i 0];
    plot_struct.k.label = {'$k_x \ (\pi/3a_{CC})$', '$k_y \ (\pi/\sqrt{3}a_{CC})$'};
    plot_struct.k.ticklabels = 0:1/(n_ticks-1):1;
    plot_struct.k.ticks = n_points*plot_struct.k.ticklabels;
end

