% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by main() to plot energy bands
% Inputs
%   E = array of energy bands in Ry ordered as (k,G)
%   config = container for axis settings

function plot_bands(E, config)
    Ef = max(E(:,4)); % set zero-point of energy
    figure();
    plot(E-Ef);
    xlabel(config.k.label);
    ylabel(config.E.label);
    axis([1, size(E,1), -config.E.lim, config.E.lim]);
    xticks(config.k.ticks);
    xticklabels(config.k.ticklabels);
    title(config.nanostructure);
    savefig(append('Figures/',config.nanostructure,'_bands.fig'));
end