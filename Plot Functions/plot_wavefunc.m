% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by main() to plot cross-sectional wavefunctions
% Inputs:
%   param = container for nanostructure parameters
%   config = container for figure/axis settings

function plot_wavefunc(param)
    unit_conv = 1e9*param.r_H;
    u = param.Bloch_comp{1};
    G = param.G_vec/unit_conv;
    R_gen = param.R_gen*unit_conv;
    V_cell = dot(cross(R_gen(:,1),R_gen(:,2)),R_gen(:,3));
    n_step = 100;
    x_max = R_gen(1,1)/2;
    x_min = -x_max;
    dx = (x_max-x_min)/n_step;
    x = x_min:dx:x_max;
    n_x = length(x);
    y_max = R_gen(2,2)/2;
    y_min = -y_max;
    dy = (y_max-y_min)/n_step;
    y = y_min:dy:y_max;
    n_y = length(y);
    z_max = R_gen(3,3)/2;
    z_min = -z_max;
    dz = (z_max-z_min)/n_step;
    z = z_min:dz:z_max;
    n_z = length(z);
    n_band = 6;
    psi = zeros(n_y,n_z,n_band);
    psi_x = zeros(n_x,1);
    figure();
    for i_band = 1:n_band
        n_band = param.n_valence-3+i_band;
        u_i = u(:,n_band);
        for i_y = 1:n_y
            y_i = y(i_y);
            for i_z = 1:n_z
                z_i = z(i_z);
                for i_x = 1:n_x
                    r = [x(i_x), y_i, z_i];
                    psi_x(i_x) = sum(u_i.*exp(1i*dot(G, r.*ones(size(G)), 2)), 1);
                end
                psi(i_y,i_z,i_band) = mean(abs(psi_x).^2)/V_cell;
            end
        end
        subplot(2,3,i_band);
        surf(y,z,abs(psi(:,:,i_band)).^2);
        xlabel('y (nm)');
        ylabel('z (nm)');
        zlabel('$|\psi|^2$');
        xlim([y_min y_max]);
        ylim([z_min z_max]);
        title(append('band ',num2str(n_band)));
        colorbar;
        view(2);
    end
    sgtitle(param.nanostructure);
    savefig(append('Figures/',param.nanostructure,'_wavefunc.fig'));
end