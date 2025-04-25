% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by select_param() to define superlattice parameters specific to zigzag nanotubes.
% Inputs:
%   N_1 = number of one primitive graphene translation for wrapping vector
%   param = container for nanostructure parameters
%   config = container for figure/axis settings
% Outputs:
%   updated param and config
%
% References:
%   [4] M. V. Fischetti and W. G. Vandenberghe. 
%       Advanced Physics of Electron Transport in Semiconductors and Nanostructures.
%       Springer (2016).

function [param, config] = zNT(N_1, param, config)

    a_CC = param.a_CC;
    a_CH = param.a_CH;
    N_a = param.N_a;
    N_x = param.N_x;
    N_z = param.N_z;

    x_hat = [1; 0; 0]; % unit vector along tube length (axial)
    r_hat = @(phi) [0; cos(phi); sin(phi)]; % radial unit vector as a function of azimuthal angle phi
    R_graphene = sqrt(3)*a_CC*[1/2 -1/2 0; sqrt(3)/2 sqrt(3)/2 0]; % primitive lattice vectors of graphene
    C = N_1*R_graphene(:,1); % wrapping vector (circumference)
    r_tube = norm(C)/(2*pi); % nanotube radius
    theta = a_CC/r_tube; % characteristic angle of carbon-carbon bond
    n_C = 4*N_1*N_a; % number of carbon atoms within supercell
    n_H = 2*N_1*(N_x>0); % number of hydrogen atoms to terminate dangling bonds
    tau = zeros(3,n_C+n_H); % start locating carbon atoms
    
    % locate carbon atoms within supercell
    idx_C = 0; % number of carbon atoms counted so far
    for l = 0:floor(N_a-1) % loop over tube length within supercell
        for c = 0:N_1-1 % loop over tube circumference
            tau(:,idx_C+1) = r_tube*r_hat((c+1/2)*sqrt(3)*theta) + l*3*a_CC*x_hat;
            tau(:,idx_C+2) = r_tube*r_hat(c*sqrt(3)*theta) + (l*3+1/2)*a_CC*x_hat;
            tau(:,idx_C+3) = r_tube*r_hat(c*sqrt(3)*theta) + (l*3+3/2)*a_CC*x_hat;
            tau(:,idx_C+4) = r_tube*r_hat((c+1/2)*sqrt(3)*theta) + (l*3+2)*a_CC*x_hat;
            idx_C = idx_C + 4;
        end
    end
    if mod(N_a,1) > 0 % half-integer N_a; add half-motif along circumference
        for c = 0:N_1-1
            tau(:,idx_C+1) = tau(:,4*c+1) + floor(N_a)*3*a_CC*x_hat;
            tau(:,idx_C+2) = tau(:,4*c+2) + floor(N_a)*3*a_CC*x_hat;
            idx_C = idx_C + 2;
        end
    end

    if n_H > 0 % finite tubes; locate hydrogen atoms along edges within supercell
        idx_H = 0; % number of hydrogen atoms counted so far
        for c = 0:N_1-1
            tau(:,n_C+idx_H+1) = tau(:,4*c+1) - a_CH*x_hat;
            tau(:,n_C+idx_H+2) = (mod(N_a,1)>0)*(tau(:,4*c+2) + (floor(N_a)*3*a_CC+a_CH)*x_hat) + (mod(N_a,1)==0)*(tau(:,4*c+4) + (floor(N_a-1)*3*a_CC+a_CH)*x_hat);
            idx_H = idx_H + 2;
        end
    end

    param.tau = tau - ones(size(tau)).*mean(tau,2); % shift coordinate origin to center of nanostructure

    % atomic radii
    param.r_atom = ones(1,size(tau,2));
    param.r_atom(1:n_C) = param.r_C;

    % generators of superlattice vectors
    tau_range = range(tau,2);
    param.R_gen = [(N_x>0)*(tau_range(1)+N_x*3*a_CC)+(N_x==0)*3*a_CC 0 0; 0 tau_range(2)+N_z*sqrt(3)*a_CC 0; 0 0 tau_range(3)+N_z*sqrt(3)*a_CC];

    % real-space path along tube axis
    param.R = (0:1/(config.n_points-1):1)'*(param.R_gen(:,1))';
    config.R.label = {'$x \ (3a_{CC})$'};
    config.R.ticklabels = 0:1/(config.n_ticks-1):1;
    config.R.ticks = config.n_points*config.R.ticklabels;

    % reciprocal-space path along tube axis
    k_max = pi/vecnorm(param.R_gen(:,1));
    param.dk = k_max/(config.n_points-1);
    param.k = (0:param.dk:k_max)'*[1 0 0];
    config.k.label = {'$k_x \ (\pi/3a_{CC})$'};
    config.k.ticklabels = 0:1/(config.n_ticks-1):1;
    config.k.ticks = config.n_points*config.k.ticklabels;
end