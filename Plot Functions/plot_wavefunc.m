% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by nano_C to plot energy bands
% Inputs:
%   psi = cell of wavefunctions ordered by increasing eigenenergy
%   structure = name of selected nanostructure
%   plot_struct = container for axis settings

function plot_wavefunc(psi, structure, plot_struct)
    R = (1:size(psi,1))';
    k = (1:size(psi,2))';
    interp = 'latex';
    figure();
    hold on
    for n = 1:size(psi,3)
        mesh(R, k, abs(psi(:,:,n)).^2);
    end
    hold off
    xlabel(plot_struct.R.label,'Interpreter',interp);
    ylabel(plot_struct.k.label,'Interpreter',interp);
    zlabel('$|\psi|^2$','Interpreter',interp);
    set(gca,'TickLabelInterpreter',interp);
    xticks(plot_struct.R.ticks);
    xticklabels(plot_struct.R.ticklabels);
    yticks(plot_struct.k.ticks);
    yticklabels(plot_struct.k.ticklabels);
    title(structure);
    savefig(append('Figures/',structure,'_wavefunc.fig'));
end