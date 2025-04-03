% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by nano_C to compute energy bands.
% Inputs:
%   k = reciprocal-space points at which eigenenergies are computed
%   R_gen = generators of superlattice vectors
%   R = lattice translation over which to evaluate wavefunctions
%   tau = atomic positions within the supercell
%   r_atom = atomic radii ordered as tau
%   a_d = diamond lattice constant in Bohr radii
%   E_cut = cutoff energy in Ry
%   max_G = maximum manhattan distance of reciprocal lattice vectors
% Output:
%   E = array of energy bands in Ry ordered as (k,G)
%   psi = cell of wavefunctions ordered by increasing eigenenergy
%
% References:
%   [1] W. G. Vandenberghe. 
%       bulk_pseudo_william.m.
%   [2] Y. Kurokawa, S. Nomura, T. Takemori, and Y. Aoyagi.
%       "Large-scale calculation of optical dielectric functions of diamond nanocrystallites."
%       Phys. Rev. B 61 12616 (2000).
%   [3] M. V. Fischetti and W. G. Vandenberghe. 
%       Advanced Physics of Electron Transport in Semiconductors and Nanostructures.
%       Springer (2016).

function [E, psi] = bands(k, R_gen, R, tau, r_atom, a_d, E_cut, max_G)

    % empirical pseudopotentials [2]
    b_C = [1.781 1.424 0.354 0.938]; % empirical parameters for carbon
    U_C = @(q_in) (a_d/2)^3*(b_C(1)*(b_C(3)*q_in.^2-b_C(2))./(exp(b_C(3)*q_in.^2-b_C(4))+1)); % pseudopotential function for carbon
    b_H = [-0.397 0.0275 0.1745 -0.0531 0.0811 -1.086 2.71 -2.86]; % empirical parameters for hydrogen
    U_H = @(q_in) (a_d/2)^3*((q_in<=2).*(b_H(1) + b_H(2)*q_in + b_H(3)*q_in.^2 + b_H(4)*q_in.^3) + (q_in>2).*(b_H(5)./(q_in+(q_in==0)) + b_H(6)./(q_in.^2+(q_in==0)) + b_H(7)./(q_in.^3+(q_in==0)) + b_H(8)./(q_in.^4+(q_in==0)))); % pseudopotential function for hydrogen

    n_k = size(k,1); % number of k points
    n_atoms = size(tau,2); % number of atoms within supercell
    V_cell = dot(cross(R_gen(:,1),R_gen(:,2)),R_gen(:,3)); % supercell volume

    % generate reciprocal lattice vectors [1]
    G_gen = (2*pi/V_cell)*[cross(R_gen(:,2),R_gen(:,3)), cross(R_gen(:,3),R_gen(:,1)), cross(R_gen(:,1),R_gen(:,2))];
    combo_gen = (-max_G:max_G)';
    G = [ kron(kron(combo_gen,combo_gen.^0),combo_gen.^0) kron(kron(combo_gen.^0,combo_gen),combo_gen.^0) kron(kron(combo_gen.^0,combo_gen.^0),combo_gen) ]*G_gen';
    G = G(G.^2*ones(3,1)<E_cut,:); % retain all G vectors within the cut-off
    n_G = size(G,1);
    GmG2 = kron(ones(1,n_G), G) - kron(ones(n_G,1), reshape(G',1,3*n_G)); % compute differences between G and G' (size 3*nG^2) [1]
    
    % off-diagonal (potential) elements of Hamiltonian
    U = zeros(n_G);
    disp('Computing pseudopotentials...')
    tic
    for i_atom = 1:n_atoms
        if r_atom(i_atom) == 1 % hydrogen
            U = U + (1/V_cell)*exp(1i*(GmG2)*kron(eye(n_G),tau(:,i_atom))).*U_H(sqrt(abs(GmG2).^2*kron(eye(n_G),ones(3,1))));
        else % carbon
            U = U + (1/V_cell)*exp(1i*(GmG2)*kron(eye(n_G),tau(:,i_atom))).*U_C(sqrt(abs(GmG2).^2*kron(eye(n_G),ones(3,1))));
        end
    end
    toc
    E = zeros(n_k,n_G);
    u = cell(n_k,1);
    disp('Computing eigenvalues...')
    tic
    for i_k = 1:n_k
        T = (spdiags(sum(abs(G+ones(n_G,1)*k(i_k,:)).^2, 2), 0, n_G, n_G)); % diagonal (kinetic) elements of Hamiltonian
        [u_mat, E_mat] = eig(T+U); % eigenvalues of Hamiltonian
        E(i_k,:) = sort(real(ones(1,n_G)*E_mat));
        u{i_k} = sort(u_mat, 1);
    end
    toc
    n_psi = 1; % number of wavefunctions to compute, starting from the lowest band
    n_R = size(R,1); % number of real-space points
    psi = zeros(n_R,n_k,n_psi);
    disp('Computing wavefunctions...')
    tic
    for i_psi = 1:n_psi
        for i_k = 1:n_k
            for i_R = 1:n_R
                psi(i_R,i_k,i_psi) = (1/V_cell)*sum(u{i_k}(:,i_psi).*exp(1i*dot(k(i_k,:).*ones(size(G))+G, R(i_R,:).*ones(size(G)), 2)), 1);
            end
        end
    end
    toc
    disp('Finished.');
end