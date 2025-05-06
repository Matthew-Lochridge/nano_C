% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by main() to plot energy bands
% Inputs
%   param = container for nanostructure parameters
%   config = container for figure/axis settings

function plot_bands(param, config)
    r_H = param.r_H;
    Ry = param.Ry;
    E_band = param.E*Ry;
    grad_E_band = param.grad_E*Ry*r_H*100;
    grad_E_band(grad_E_band==0) = realmin;
    E_max = config.E.lim;
    E_min = -E_max;
    E = E_min:5e-3*(E_max-E_min)/(size(E_band,2)-1):E_max;
    DOS = zeros(size(E));
    G_ball = zeros(size(E));
    abs_tol = 1e-1;
    parfor i_E = 1:length(E)
        for i_band = 1:size(E_band,2)
            G_ball(i_E) = G_ball(i_E) + ((sum(E_band(:,i_band)<=E(i_E),'all')>0) && (sum(E_band(:,i_band)>=E(i_E),'all')>0));
        end
        DOS(i_E) = sum(isapprox(E(i_E)*ones(size(E_band)), E_band, AbsoluteTolerance=abs_tol)./abs(grad_E_band),'all')/pi;
    end
    figure();
    subplot(1,3,1);
    plot(E_band,'-k');
    xlabel(config.k.label);
    ylabel('$E$ (eV)');
    axis([1, size(E_band,1), E_min, E_max]);
    xticks(config.k.ticks);
    xticklabels(config.k.ticklabels);
    subplot(1,3,2);
    stairs(1e-8*DOS,E,'-k');
    xlabel('DOS ($10^{8}$ eV$^{-1}$ cm$^{-1}$)');
    ylim([E_min,E_max]);
    subplot(1,3,3);
    stairs(G_ball,E,'-k');
    xlabel('G ($2e^2/h$)');
    ylim([E_min,E_max]);
    sgtitle(param.nanostructure);
    savefig(append('Figures/',param.nanostructure,'_bands.fig'));
end