% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by main() to compute energy bands.
% Inputs:
%   param = container for nanostructure parameters
%   U_C = empirical pseudopotential function handle for carbon atoms
%   U_H = empirical pseudopotential function handle for hydrogen atoms
% Output:
%   E = array of energy bands in Ry ordered as (k,G)
%   psi = cell of wavefunctions ordered by increasing eigenenergy
%
% References:
%   [1] W. G. Vandenberghe. 
%       bulk_pseudo_william.m.
%   [4] M. V. Fischetti and W. G. Vandenberghe. 
%       Advanced Physics of Electron Transport in Semiconductors and Nanostructures.
%       Springer (2016).

function [E, psi] = bands(param, U_C, U_H)

    k = param.k;
    R = param.R;
    R_gen = param.R_gen;
    tau = param.tau;
    r_atom = param.r_atom;

    n_k = size(k,1); % number of reciprocal-space points
    n_R = size(R,1); % number of real-space points
    n_atoms = size(tau,2); % number of atoms within supercell
    V_cell = dot(cross(R_gen(:,1),R_gen(:,2)),R_gen(:,3)); % supercell volume

    % generate reciprocal lattice vectors [1]
    G_gen = (2*pi/V_cell)*[cross(R_gen(:,2),R_gen(:,3)), cross(R_gen(:,3),R_gen(:,1)), cross(R_gen(:,1),R_gen(:,2))];
    combo_gen = (-param.max_G:param.max_G)';
    G = [ kron(kron(combo_gen,combo_gen.^0),combo_gen.^0) kron(kron(combo_gen.^0,combo_gen),combo_gen.^0) kron(kron(combo_gen.^0,combo_gen.^0),combo_gen) ]*G_gen';
    G = G(G.^2*ones(3,1)<param.E_cut,:); % retain all G vectors within the cut-off
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
    parfor i_k = 1:n_k
        T = (spdiags(sum(abs(G+ones(n_G,1)*k(i_k,:)).^2, 2), 0, n_G, n_G)); % diagonal (kinetic) elements of Hamiltonian
        [u_mat, E_mat] = eig(T+U); % eigenvalues of Hamiltonian
        E(i_k,:) = sort(real(ones(1,n_G)*E_mat));
        u{i_k} = sort(u_mat, 1);
    end
    toc
    n_psi = 1; % number of wavefunctions to compute, starting from the lowest band
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