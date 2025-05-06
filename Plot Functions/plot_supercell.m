% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by main() to plot supercells for selected nanostructures.
% Inputs:
%   param = container for nanostructure parameters
%   config = container for figure/axis settings

function plot_supercell(param)
    config = supercell_config(param);
    r_H = param.r_H*10^9;
    tau = r_H*param.tau;
    r_atom = param.r_atom;
    figure();
    hold on;
    plot3(tau(1,r_atom>1), tau(2,r_atom>1), tau(3,r_atom>1), config.atom.C);
    plot3(tau(1,r_atom==1), tau(2,r_atom==1), tau(3,r_atom==1), config.atom.H);
    for i = 1:size(tau,2)-1
        for j = i+1:size(tau,2)
            separation = tau(:,j) - tau(:,i);
            if strcmpi(param.nanostructure,'ppp')
                if isapprox(vecnorm(separation), r_H*param.a_CC_PPP_hi, AbsoluteTolerance=config.abs_tol) || isapprox(vecnorm(separation), r_H*param.a_CC_PPP_lo, AbsoluteTolerance=config.abs_tol)
                    plot3([tau(1,i), tau(1,j)], [tau(2,i), tau(2,j)], [tau(3,i), tau(3,j)], config.bond.C);
                elseif isapprox(vecnorm(separation), r_H*param.a_CH_PPP, AbsoluteTolerance=config.abs_tol)
                    plot3([tau(1,i), tau(1,j)], [tau(2,i), tau(2,j)], [tau(3,i), tau(3,j)], config.bond.H);
                end
            else
                if isapprox(vecnorm(separation), r_H*param.a_CC_graphene, AbsoluteTolerance=config.abs_tol)
                    plot3([tau(1,i), tau(1,j)], [tau(2,i), tau(2,j)], [tau(3,i), tau(3,j)], config.bond.C);
                elseif isapprox(vecnorm(separation), r_H*param.a_CH_graphene, AbsoluteTolerance=config.abs_tol)
                    plot3([tau(1,i), tau(1,j)], [tau(2,i), tau(2,j)], [tau(3,i), tau(3,j)], config.bond.H);
                end
            end
        end
    end
    xlabel(config.xlabel);
    ylabel(config.ylabel);
    zlabel(config.zlabel);
    axis(config.limits);
    title(param.nanostructure);
    legend(config.legend);
    hold off;
end

function config = supercell_config(param)
    config = struct();
    config.limits = 10^9*param.r_H*max(vecnorm(param.R_gen))*[-1 1 -1 1 -1 1]/2; % kron(vecnorm(10^9*r_H*param.R_gen)',[-1 1])/2;
    config.abs_tol = 1e-2;
    config.title = append(param.nanostructure,' supercell');
    config.xlabel = '$x$ (nm)';
    config.ylabel = '$y$ (nm)';
    config.zlabel = '$z$ (nm)';
    config.atom.C = '.k';
    config.atom.H = 'ok';
    config.bond.C = '-k';
    config.bond.H = '--k';
    config.legend = {'C'};
    if sum(param.r_atom==1) > 0 % hydrogen included
        config.legend{end+1} = 'H';
    end
end