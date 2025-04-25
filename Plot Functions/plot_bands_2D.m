% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by main() to plot energy bands in 2D reciprocal-space
% Inputs
%   param = container for nanostructure parameters
%   config = container for figure/axis settings

function plot_bands_2D(param, config)
    k = param.k;
    E_band = param.E_band*param.Ry;
    if strcmpi(config.nanostructure,'graphene') 
        [k,E_band] = reflect(k,E_band,config.vertex,config.sym_proj); % reflect IBZ into full hexagonal BZ
    end
    figure();
    subplot(1,2,1);
    plot3(k(:,1),k(:,2),E_band);
    hold off
    xlabel(config.k.label{1});
    ylabel(config.k.label{2});
    zlabel(config.E.label);
    zlim([-config.E.lim, config.E.lim]);
    text(config.k.sym_points(:,1),config.k.sym_points(:,2),config.k.sym_points(:,3),config.k.text);
    subplot(1,2,2);
    stairs(12*param.DOS/((100*param.r_H)^2*param.Ry),param.E_DOS*param.Ry);
    xlabel('DOS (eV$^{-1}$ cm$^{-2}$)');
    ylim([-config.E.lim,config.E.lim]);
    sgtitle(param.nanostructure);
    savefig(append('Figures/',param.nanostructure,'_BZ.fig'));
end

function [k, E] = reflect(k, E, v, u)
    k = [k; k + 2*point_to_line(k,zeros(1,3),v(1,:))*u(1,:)];
    k_ref_2 = k + 2*point_to_line(k,zeros(1,3),v(2,:))*u(2,:);
    k_ref_3 = k + 2*point_to_line(k,zeros(1,3),v(3,:))*u(3,:);
    k = [k; k_ref_2; k_ref_3];
    k = [k; k + 2*point_to_line(k,v(4,:),-v(4,:))*u(4,:)];
    E = [E;E;E;E;E;E;E;E;E;E;E;E];
end

function d = point_to_line(pt, v1, v2)
% pt should be nx3
% v1 and v2 are vertices on the line (each 1x3)
% d is a nx1 vector with the orthogonal distances
    a = v1 - v2;
    b = pt - ones(size(pt,1),1)*v2;
    d = vecnorm(cross(ones(size(pt,1),1)*a,b),2,2) / vecnorm(a);
end