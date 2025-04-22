% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by main() to plot energy bands
% Inputs
%   param = container for nanostructure parameters
%   config = container for figure/axis settings

function plot_bands(param, config)
    [x_DOS, D_E_per_eV_per_cm] = DOS(param);
    E = -config.E.lim:config.E.lim;
    figure();
    subplot(1,2,1);
    plot(param.E_bands*param.Ry);
    xlabel(config.k.label);
    ylabel(config.E.label);
    axis([1, size(param.E_bands,1), -config.E.lim, config.E.lim]);
    xticks(config.k.ticks);
    xticklabels(config.k.ticklabels);
    subplot(1,2,2);
    stairs(sum(param.DOS(E))/(100*param.r_H*param.Ry),E);
    xlabel('DOS ($10^{8}$ eV$^{-1}$ cm$^{-1}$)');
    ylim([-config.E.lim,config.E.lim]);
    sgtitle(config.nanostructure);
    % savefig(append('Figures/',config.nanostructure,'_bands.fig'));
end

function [x_DOS, D_E_per_eV_per_cm] = DOS(param)
    hist_steps = size(param.E_bands,2);
    [D_E_per_bin_per_a,edges] = histcounts(param.E_bands,hist_steps);
    bin_lengths = diff(edges);
    x_DOS = edges(1:end-1) + bin_lengths(1)/2;
    D_E_per_eV_per_cm = 10^(-8)*D_E_per_bin_per_a*bin_lengths(1)/(sqrt(3)*param.a_CC*param.r_H*10^2*param.Ry);
end