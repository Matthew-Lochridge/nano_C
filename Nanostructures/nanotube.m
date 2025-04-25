% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by select_param() to define superlattice parameters specific to nanotubes.
% Inputs:
%   N_1, N_2 = numbers of each primitive graphene translation for wrapping vector
%   param = container for nanostructure parameters
%   config = container for figure/axis settings
% Outputs:
%   updated param and config
%
% References:
%   [5] J. W. Mintmire and C. T. White. 
%       "Electronic and structural properties of carbon nanotubes."
%       Carbon 33 7 (1995).

function [param, config] = nanotube(N_1, N_2, param, config)

    a_CC = param.a_CC;
    a_CH = param.a_CH;
    N_a = param.N_a;
    N_x = param.N_x;
    N_z = param.N_z;

    x_hat = [1; 0; 0];
    y_hat = [0; 1; 0];
    z_hat = [0; 0; 1];
    a_hat = [1/2; sqrt(3)/2; 0]; % unit vector of primitive graphene translation
    a_1 = sqrt(3)*a_CC*a_hat;
    a_2 = sqrt(3)*a_CC*(a_hat-x_hat);
    
    GCD = gcd(N_1,N_2); % greatest common divisor of N_1 and N_2
    n = (mod(N_1-N_2,3*GCD)==0)*(3*GCD) + (mod(N_1-N_2,3*GCD)~=0)*GCD;
    C = N_1*a_1 + N_2*a_2; % wrapping vector within graphene plane
    L = (sqrt(3)/n)*cross(z_hat,C); % axial period within graphene plane
    h = @(d) sign(dot(d,y_hat))*vecnorm(cross(d,C))/vecnorm(C); % primitive helical translation (axial)
    N_L = abs(round(vecnorm(L)/h(a_1+a_2))); % number of axial translations in period
    phi = @(d) 2*pi*dot(d,C)/(vecnorm(C)^2); % mapping from two-atom separation in graphene to helical rotation
    R = @(theta) [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)]; % matrix of rotations about tube axis by angle theta
    S = @(tau, d) R(phi(d))*tau + h(d)*x_hat; % screw operation mapping translation d relative to tau from planar graphene onto nanotube
    
    % locate atoms within supercell
    n_C = 2*GCD*N_L*N_a; % number of carbon atoms within supercell
    n_H = 0; % 6*GCD*(N_x>0); % number of hydrogen atoms to terminate edge bonds
    idx_C = 0; % number of carbon atoms counted so far
    tau = zeros(3,n_C+n_H); % start counting carbon atoms
    tau(:,1) = (vecnorm(C)/(2*pi))*y_hat;
    for a = 0:floor(N_a-1) % loop over whole cells
        for l = 0:N_L-1 % loop over length of primitive helical motifs within cell
            for r = 0:GCD-1 % loop over primitive helical motif
                tau(:,idx_C+1) = R(2*pi*r/GCD)*S(tau(:,1), l*(a_1+a_2)) + a*vecnorm(L)*x_hat;
                tau(:,idx_C+2) = S(tau(:,idx_C+1), (a_1+a_2)/3);
                idx_C = idx_C + 2;
            end
        end
    end
    if mod(N_a,1) > 0 % non-integer N_a; add a half-cell along the tube
        for t = 1:idx_C/(2*floor(N_a))
            tau(:,idx_C+1) = tau(:,t) + floor(N_a)*vecnorm(L)*x_hat;
            idx_C = idx_C + 1;
        end
    end
    if n_H > 0 % start counting hydrogen atoms
        idx_H = 0; % number of hydrogen atoms counted so far
        d_H = @(d) a_CH*d/vecnorm(d); % separation between H and nearest C atoms
        for r = 0:GCD-1
            tau(:,n_C+idx_H+1) = S(tau(:,2*r+1), d_H(a_1-2*(a_1+a_2)/3));
            tau(:,n_C+idx_H+2) = S(tau(:,n_C+idx_H+1), (1+a_CH/a_CC)*(a_1+a_2)/3);
            tau(:,n_C+idx_H+3) = S(tau(:,n_C+idx_H+1), (a_1+a_2));
            %tau(:,n_C+idx_H+3) = S( (mod(N_a,1)==0)*(tau(:,2*GCD*floor(N_a)*(N_L-1)+2*r+2)+floor(N_a-1)*N_L*vecnorm(L)*x_hat) + (mod(N_a,1)>0)*(tau(:,2*GCD+2*r+2)+floor(N_a)*N_L*vecnorm(L)*x_hat), d_H(a_1,1));
            idx_H = idx_H + 3; % (10,5) only
        end
    end

    param.tau = tau - ones(size(tau)).*mean(tau,2); % shift coordinate origin to center of nanostructure

    % atomic radii
    param.r_atom = ones(1,size(tau,2));
    param.r_atom(1:n_C) = param.r_C;

    % generators of superlattice vectors
    tau_range = range(tau,2);
    param.R_gen = [(N_x>0)*(tau_range(1)+N_x*a_CC)+(N_x==0)*N_a*vecnorm(L) 0 0; 0 tau_range(2)+N_z*a_CC 0; 0 0 tau_range(3)+N_z*a_CC];

    % real-space path along tube axis
    param.R = (0:1/(config.n_points-1):1)'*(param.R_gen(:,1))';
    config.R.label = {'$x \ (a_x)$'};
    config.R.ticklabels = 0:1/(config.n_ticks-1):1;
    config.R.ticks = config.n_points*config.R.ticklabels;

    % reciprocal-space path along tube axis
    k_max = pi/vecnorm(param.R_gen(:,1));
    param.dk = k_max/(config.n_points-1);
    param.k = (0:param.dk:k_max)'*[1 0 0];
    config.k.label = {'$k_x \ (\pi/a_x)$'};
    config.k.ticklabels = 0:1/(config.n_ticks-1):1;
    config.k.ticks = config.n_points*config.k.ticklabels;
end