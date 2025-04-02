% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by nano_C to plot supercells for selected nanostructures.
% Inputs:
%   structure = name of selected nanostructure
%   a_CC = carbon-carbon bond length in nm
%   a_CH = carbon-hydrogen bond length in nm
%   tau = atom positions within the supercell in nm
%   r_atom = atomic radii in Bohr radii ordered as tau
%   N_x = separation between finite ribbons or tubes in nm (axial)
%   N_y = separation between finite ribbons or tubes in nm (transverse, in-plane for ribbons)
%   N_z = separation between ribbons or tubes in nm (transverse, out-of-plane for ribbons)

function plot_supercell(structure, a_CC, a_CH, tau, r_atom, N_x, N_y, N_z)
    supercell_range = range(tau,2) + [N_x; N_y; N_z];
    limits = kron(supercell_range,[-1 1])/2;
    figure(); xlabel('x (nm)'); ylabel('y (nm)'); zlabel('z (nm)'); hold on;
    plot3(tau(1,r_atom>1), tau(2,r_atom>1), tau(3,r_atom>1), '.k');
    plot3(tau(1,r_atom==1), tau(2,r_atom==1), tau(3,r_atom==1), 'ok');
    axis([limits(1,:) limits(2,:) limits(3,:)]);
    for k1 = 1:size(tau,2)-1
        for k2 = k1+1:size(tau,2)
            sep = tau(:,k2) - tau(:,k1);
            if isapprox(norm(sep),a_CC,AbsoluteTolerance=1e-2)
                plot3([tau(1,k1), tau(1,k2)], [tau(2,k1), tau(2,k2)], [tau(3,k1), tau(3,k2)], '-k');
            elseif isapprox(norm(sep),a_CH,AbsoluteTolerance=1e-2)
                plot3([tau(1,k1), tau(1,k2)], [tau(2,k1), tau(2,k2)], [tau(3,k1), tau(3,k2)], '--k');
            end
        end
    end
    hold off; title(append(structure,' supercell')); legend('C','H');
end