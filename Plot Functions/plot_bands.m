% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by nano_C to plot energy bands
% Inputs
%   E = array of energy bands in Ry ordered as (k,G)
%   E_max = maximum energy shown in band plot in eV
%   structure = name of selected nanostructure
%   plot_struct = container for axis settings

function plot_bands(E, E_max, structure, plot_struct)
    interp = 'latex';
    figure();
    Ef = max(E(:,4)); % set zero-point of energy
    plot(E-Ef);
    xlabel(plot_struct.k.label,'Interpreter',interp);
    ylabel('$E$ (eV)','Interpreter',interp);
    axis([1 size(E,1) -E_max E_max]);
    set(gca,'TickLabelInterpreter',interp);
    xticks(plot_struct.k.ticks);
    xticklabels(plot_struct.k.ticklabels);
    title(structure);
    savefig(append('Figures/',structure,'_bands.fig'));
end