% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by main() to plot energy bands
% Inputs:
%   psi = cell of wavefunctions ordered by increasing eigenenergy
%   config = container for figure/axis settings

function plot_wavefunc(psi, config)
    R = (1:size(psi,1))';
    k = (1:size(psi,2))';
    figure();
    hold on
    for n = 1:size(psi,3)
        mesh(R, k, abs(psi(:,:,n)).^2);
    end
    hold off
    xlabel(config.R.label);
    ylabel(config.k.label);
    zlabel(config.psi2.label);
    xticks(config.R.ticks);
    xticklabels(config.R.ticklabels);
    yticks(config.k.ticks);
    yticklabels(config.k.ticklabels);
    title(config.nanostructure);
    savefig(append('Figures/',config.nanostructure,'_wavefunc.fig'));
end