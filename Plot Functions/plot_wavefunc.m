% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by main() to plot energy bands
% Inputs:
%   param = container for nanostructure parameters
%   config = container for figure/axis settings

function plot_wavefunc(param, config)
    psi = param.psi;
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
    xlim([R(1) R(end)]);
    ylim([k(1) k(end)]);
    title(param.nanostructure);
    savefig(append('Figures/',param.nanostructure,'_wavefunc.fig'));
end