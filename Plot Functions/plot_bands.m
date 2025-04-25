% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by main() to plot energy bands
% Inputs
%   param = container for nanostructure parameters
%   config = container for figure/axis settings

function plot_bands(param, config)
    E_band = param.E_band*param.Ry;
    figure();
    subplot(1,2,1);
    plot(E_band);
    xlabel(config.k.label);
    ylabel(config.E.label);
    axis([1, size(E_band,1), -config.E.lim, config.E.lim]);
    xticks(config.k.ticks);
    xticklabels(config.k.ticklabels);
    subplot(1,2,2);
    stairs(1e-8*param.DOS/((100*param.r_H)*param.Ry),param.E_DOS);
    xlabel('DOS ($10^{8}$ eV$^{-1}$ cm$^{-1}$)');
    ylim([-config.E.lim,config.E.lim]);
    sgtitle(param.nanostructure);
    savefig(append('Figures/',param.nanostructure,'_bands.fig'));
end