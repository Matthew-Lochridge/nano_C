% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by main() to plot supercells for selected nanostructures.
% Inputs:
%   param = container for nanostructure parameters
%   config = container for figure/axis settings

function plot_supercell(param, config)
    r_H = param.r_H*10^9;
    tau = r_H*param.tau;
    r_atom = param.r_atom;
    limits = max(vecnorm(r_H*param.R_gen))*[-1 1 -1 1 -1 1]/2; % kron(vecnorm(r_H*param.R_gen)',[-1 1])/2;
    figure();
    hold on;
    plot3(tau(1,r_atom>1), tau(2,r_atom>1), tau(3,r_atom>1), config.atom.C);
    plot3(tau(1,r_atom==1), tau(2,r_atom==1), tau(3,r_atom==1), config.atom.H);
    for i = 1:size(tau,2)-1
        for j = i+1:size(tau,2)
            separation = tau(:,j) - tau(:,i);
            if isapprox(vecnorm(separation), r_H*param.a_CC, AbsoluteTolerance=config.abs_tol)
                plot3([tau(1,i), tau(1,j)], [tau(2,i), tau(2,j)], [tau(3,i), tau(3,j)], config.bond.C);
            elseif isapprox(vecnorm(separation), r_H*param.a_CH, AbsoluteTolerance=config.abs_tol)
                plot3([tau(1,i), tau(1,j)], [tau(2,i), tau(2,j)], [tau(3,i), tau(3,j)], config.bond.H);
            end
        end
    end
    hold off;
    xlabel(config.xlabel);
    ylabel(config.ylabel);
    zlabel(config.zlabel);
    axis(limits);
    legend(config.legend);
    title(config.title);
end