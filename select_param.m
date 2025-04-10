% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
% 
% Function called by main() to select specified nanostructure parameters to generate.
% Inputs:
%   param = container for nanostructure parameters
%   config = container for figure/axis settings
% Outputs:
%   updated param and config

function [param, config] = select_param(param, config)
    addpath('Nanostructures\'); % enable use of functions to generate specified nanostructure parameters

    [param, allotrope] = parse_input(param, config.nanostructure); % decompose input from main() into parameters and allotrope

    % define k, tau, r_atom, R_gen, and plot settings based on type and parameters
    if strcmpi(allotrope,'graphene')
        disp('Graphene selected (default).')
        [param, config] = graphene(param, config);

    elseif isfield(param,'N_dw') && isfield(param,'N_dl') % nanoribbons
        disp(append(config.nanostructure,' selected.'))
        if (param.N_dl == 2) && strcmpi(allotrope,'agnr')
            [param, config] = aNR(param, config);
        elseif (param.N_dl == 2) && strcmpi(allotrope,'zgnr')
            [param, config] = zNR(param, config);
        elseif strcmpi(allotrope,'gnr')
            [param, config] = nanoribbon(param, config);
        end

    elseif isfield(param,'N_1') && isfield(param,'N_2') % nanotubes
        disp(append(config.nanostructure,' selected.'))
        k = [];
        tau = [];
        r_atom = [];
        R_gen = [];
        R = [];
        N_1 = param.N_1;
        N_2 = param.N_2;
        circ = zeros(size(N_1));
        circ_func = @(n,m) sqrt(n^2+m^2+n*m); % circumference of (n,m)-CNT in sqrt(3)*a_CC
        for t = 1:length(N_1) % loop over coaxial tubes
            if N_1(t) == 0 || N_2(t) == 0 % zigzag
                [param_t, config_t] = zNT(max(abs([N_1(t) N_2(t)])), param, config);
            elseif N_2(t) == N_1(t) % armchair
                [param_t, config_t] = aNT(N_1(t), param, config);
            else % chiral
                [param_t, config_t] = nanotube(N_1(t), N_2(t), param, config);
            end
            % append atoms to supercell for each tube
            tau = [tau, param_t.tau];
            r_atom = [r_atom, param_t.r_atom];
            % use only k, R, R_gen, and config for the largest tube
            circ(t) = circ_func(N_1(t),N_2(t));
            if circ(t) == max(circ)
                k = param_t.k;
                R = param_t.R;
                R_gen = param_t.R_gen;
                config = config_t;
            end
        end
        param.k = k;
        param.tau = tau;
        param.r_atom = r_atom;
        param.R_gen = R_gen;
        param.R = R;
    end
end

