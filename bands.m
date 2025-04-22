% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
%
% Function called by main() to compute energy bands.
% Inputs:
%   param = container for nanostructure parameters
%   U_C = empirical pseudopotential function handle for carbon atoms
%   U_H = empirical pseudopotential function handle for hydrogen atoms
% Output:
%   updated param
%
% References:
%   [1] W. G. Vandenberghe, bulk_pseudo_william.m (UT Dallas)
%   [4] M. V. Fischetti and W. G. Vandenberghe. 
%       Advanced Physics of Electron Transport in Semiconductors and Nanostructures.
%       Springer (2016).
%   [7] P. Hadley, cnt.m (TU Graz)

function param = bands(param, U_C, U_H)

    k = param.k;
    dk = param.dk;
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
    param.G = G;
    param.GmG2 = GmG2;

    % compute energy bands
    U = zeros(n_G); % off-diagonal (potential) elements of Hamiltonian
    disp('Computing pseudopotentials...');
    tic
    parfor i_atom = 1:n_atoms
        if r_atom(i_atom) == 1 % hydrogen
            U = U + (1/V_cell)*exp(1i*(GmG2)*kron(eye(n_G),tau(:,i_atom))).*U_H(sqrt(abs(GmG2).^2*kron(eye(n_G),ones(3,1))));
        else % carbon
            U = U + (1/V_cell)*exp(1i*(GmG2)*kron(eye(n_G),tau(:,i_atom))).*U_C(sqrt(abs(GmG2).^2*kron(eye(n_G),ones(3,1))));
        end
    end
    toc
    E = zeros(n_k,n_G);
    u = cell(n_k,1);
    grad_E = cell(n_k,1);
    dos = @(E_var) 0;
    disp('Computing energy bands and DoS...');
    tic
    parfor i_k = 1:n_k
        k_i = k(i_k,:);
        T = (spdiags(sum(abs(G+ones(n_G,1)*k_i).^2, 2), 0, n_G, n_G)); % diagonal (kinetic) elements of Hamiltonian
        [u_mat, E_mat] = eig(T+U); % eigenvalues of Hamiltonian
        [E(i_k,:),sort_idx] = sort(real(ones(1,n_G)*E_mat));
        u{i_k} = u_mat(:,sort_idx);
        grad_T = (spdiags(sum(2*(G(:,1)+ones(n_G,1)*k_i(1)+1i*(G(:,2)+ones(n_G,1)*k_i(2))), 2), 0, n_G, n_G));
        i_dE = [real(eig(real(grad_T))), real(eig(imag(grad_T))), zeros(n_G,1)];
        grad_E{i_k} = i_dE(sort_idx,:);
    end
    n_C = sum(r_atom>1); % number of carbon atoms
    n_H = sum(r_atom==1); % number of hydrogen atoms
    param.n_valence = (4*n_C+n_H)/2; % number of valence bands
    param.E_bands = E - max(E(:,param.n_valence)); % rescale energy bands
    param.u = u;
    param.grad_E = grad_E;
    abs_dE = vecnorm(grad_E,2,2);
    alpha = acos(grad_E(:,1)./abs_dE);
    w0 = (dk/2)*(cos(alpha)-sin(alpha));
    w1 = w0 + dk*sin(alpha);
    w = @(E_var) (E_var*ones(size(E_center))-E_center)./abs_dE;
    L = @(E_var) (w(E_var)<=w0).*(dk./cos(alpha)) + (w(E_var)>=w0 && w(E_var)<=w1).*(w1-w(E_var))./(cos(alpha).*sin(alpha));
    dos = @(E_var) dos(E_var) + (1/(2*pi^2))*sum(L(E_var)./abs_dE);
    toc

    % compute wavefunctions
    n_psi = 1; % number of wavefunctions to compute, starting from the lowest band
    psi = zeros(n_R,n_k,n_psi);
    disp('Computing wavefunctions...');
    tic
    for i_psi = 1:n_psi
        parfor i_k = 1:n_k
            for i_R = 1:n_R
                psi(i_R,i_k,i_psi) = (1/V_cell)*sum(u{i_k}(:,i_psi).*exp(1i*dot(k(i_k,:).*ones(size(G))+G, R(i_R,:).*ones(size(G)), 2)), 1);
            end
        end
    end
    param.psi = psi;
    toc
end