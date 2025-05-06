% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by select_param() to define superlattice parameters specific to armchair nanoribbons (infinite length).
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

function [param, config] = aNR(param, config)

    a_CC = param.a_CC_graphene;
    a_CH = param.a_CH_graphene;

    x_hat = [1; 0; 0]; % unit vector along ribbon length
    y_hat = [0; 1; 0]; % unit vector along ribbon width
    a_hat = [1/2; sqrt(3)/2; 0]; % unit vector of primitive cell
    N_w = param.N_dw/2; % number of primitive cells across ribbon width
    n_C = 4*N_w; % number of carbon atoms within supercell
    n_H = 4; % number of hydrogen atoms to terminate dangling edge bonds
    tau = zeros(3,n_C+n_H);
    
    % locate carbon atoms within supercell
    tau(:,2) = a_CC*a_hat;
    tau(:,3) = tau(:,2) + a_CC*x_hat;
    tau(:,4) = tau(:,1) + 2*a_CC*x_hat;
    for w = 0:floor(N_w-1) % loop over whole primitive cells along ribbon width
        for n = 1:4 % loop over atoms within primitive cell
            tau(:,w*4+n) = tau(:,n) + w*sqrt(3)*a_CC*y_hat;
        end
    end
    if mod(N_w,1) ~= 0 % odd N_dw; add a half-primitive-cell along ribbon length
        tau(:,4*floor(N_w)+1) = tau(:,1) + floor(N_w)*sqrt(3)*a_CC*y_hat;
        tau(:,4*floor(N_w)+2) = tau(:,4) + floor(N_w)*sqrt(3)*a_CC*y_hat;
    end

    % locate hydrogen atoms along edges within supercell
    tau(:,n_C+1) = tau(:,1) + a_CH*(x_hat-a_hat);
    tau(:,n_C+2) = tau(:,4) + a_CH*(-a_hat);
    tau(:,n_C+3) = (mod(N_w,1)>0)*(tau(:,1)+floor(N_w)*sqrt(3)*a_CC*y_hat+a_CH*a_hat) + (mod(N_w,1)==0)*(tau(:,2)+floor(N_w-1)*sqrt(3)*a_CC*y_hat+a_CH*(a_hat-x_hat));
    tau(:,n_C+4) = (mod(N_w,1)>0)*(tau(:,4)+floor(N_w)*sqrt(3)*a_CC*y_hat+a_CH*(a_hat-x_hat)) + (mod(N_w,1)==0)*(tau(:,3)+floor(N_w-1)*sqrt(3)*a_CC*y_hat+a_CH*a_hat);

    param.tau = tau - ones(size(tau)).*mean(tau,2); % shift coordinate origin to center of nanostructure

    % atomic radii
    param.r_atom = ones(1,size(tau,2));
    param.r_atom(1:n_C) = param.r_C;

    % generators of superlattice vectors
    param.R_gen = [3*a_CC 0 0; 0 range(tau(2,:))+param.N_y*sqrt(3)*a_CC 0; 0 0 param.N_z*sqrt(3)*a_CC];

    % reciprocal-space path along ribbon axis
    k_max = pi/vecnorm(param.R_gen(:,1));
    param.dk = k_max/(config.n_points-1);
    param.k = (0:param.dk:k_max)'*[1 0 0];
    config.k.label = {'$k_x \ (\pi/3a_{CC})$'};
    config.k.ticklabels = 0:1/(config.n_ticks-1):1;
    config.k.ticks = config.n_points*config.k.ticklabels;
end