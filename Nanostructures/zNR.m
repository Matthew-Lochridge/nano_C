% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by select_param() to define superlattice parameters specific to zigzag nanoribbons (infinite length).
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

function [param, config] = zNR(param, config)

    a_CC = param.a_CC;
    a_CH = param.a_CH;

    x_hat = [1; 0; 0]; % unit vector along ribbon length
    y_hat = [0; 1; 0]; % unit vector along ribbon width
    a_hat = [sqrt(3)/2; 1/2; 0]; % unit vector of primitive cell
    N_w = param.N_dw/2; % number of primitive cells across ribbon width
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

    param.tau = tau - ones(size(tau)).*mean(tau,2); % shift coordinate origin to center of nanostructure

    % atomic radii
    param.r_atom = ones(1,size(tau,2));
    param.r_atom(1:n_C) = param.r_C;

    % generators of superlattice vectors
    param.R_gen = [sqrt(3)*a_CC 0 0; 0 range(tau(2,:))+param.N_y*3*a_CC 0; 0 0 param.N_z*3*a_CC];

    % real-space path along ribbon axis
    param.R = (0:1/(config.n_points-1):1)'*(param.R_gen(:,1))';
    config.R.label = {'$x \ (\sqrt{3}a_{CC})$'};
    config.R.ticklabels = 0:1/(config.n_ticks-1):1;
    config.R.ticks = config.n_points*config.R.ticklabels;

    % reciprocal-space path along ribbon axis
    k_max = pi/vecnorm(param.R_gen(:,1));
    dk = k_max/(config.n_points-1);
    param.k = (0:dk:k_max)'*[1 0 0];
    config.k.label = {'$k_x \ (\pi/\sqrt{3}a_{CC})$'};
    config.k.ticklabels = 0:1/(config.n_ticks-1):1;
    config.k.ticks = config.n_points*config.k.ticklabels;
end