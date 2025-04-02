function E=bulk_pseudo_william(a0, b, n_atoms, tau, E_cut)
% bulk_pseudo_william calculates and plots the bandstructure using
% pseudopotentials
%
%   E = bulk_pseudo_william(a0, b, n_atoms, tau, E_cut)
%      a0 - lattice constant in atomic units
%      b - 4 coefficients for pseudopotential
%         V_q = b(1)*( q^2-b(2) )./(b(3)*exp( b(4)*q^2 )-1);
%      n_atoms - number of atoms
%      tau - matrix made of columns indicating position of each atom
%      E_cut - cut-off energy in Rydbergs
%
%   E = bulk_pseud_william() calls bulk_pseudo_william(a0, b, n_atoms, tau,
%   E_cut) with default values
%
%   @author William Vandenberghe

    %% Constants
    ryd_m = 5.2917725e-11;
    ryd_eV = 13.6;
    
    if nargin==5 % Calculate and plot the eigenvalues

        % Define the pseudopotential function
        V_q = @(q_in) b(1)*( q_in.^2-b(2) )./(b(3)*exp( b(4)*q_in.^2 ) -1);


        %%%% Generate array of k-values for which the bandstructure has to be
        % solved
        g0=2*pi/a0;
        kmax=g0;
        dk  =.02*kmax;
        k   =[(kmax*sqrt(3)/2:-dk:0)'*[1 1 1]/sqrt(3) ; (0:dk:kmax)'*[1 0 0] ;  ...
            (0:dk:kmax/2/sqrt(2))'.^0*[kmax 0 0] + (0:dk:kmax/2/sqrt(2))'*[0 1 1]/sqrt(2) ;  ...
            (kmax*3*sqrt(2)/4:-dk:0)'*[0 1 1]/sqrt(2)];
        nk=size(k,1); % Number of k-points
        %%%%

        %%%% Generate all relevant reciprocal vectors
        % Generators of the group of reciprocal vectors for FCC
        Ggen=g0*[ 1 1 -1 ; 1 -1 1 ; -1 1 1 ];

        % Generate all linear combinations of the generators within manhattan
        % distance 2
        comb=(-2:2)';
        G=[ kron(kron(comb,comb.^0),comb.^0) kron(kron(comb.^0,comb),comb.^0) kron(kron(comb.^0,comb.^0),comb) ]*Ggen;

        % Retain all reciprocal vectors within the cut-off
        G = G(G.^2*ones(3,1)<E_cut,:);
        % Number of G-vectors
        nG = size(G,1);
        %%%%


        % Generate the potential part of the Hamiltonian
        GmG2=(kron(ones(1,nG), G)-kron(ones(nG,1), reshape(G',1,3*nG))); % Matrix containing differences between G and G' (size 3*nG^2)
        A2=zeros(nG); % Initialize potential part of the Hamiltionian
        for i_atom=1:n_atoms
            A2=A2+1/n_atoms*exp(1i*(GmG2)*kron(eye(nG),tau(:,i_atom))).*V_q(sqrt(abs(GmG2).^2*kron(eye(nG),ones(3,1))));
        end


        E=zeros(nk,nG);  % Initialize the eigenvalues
        for i_k=1:nk
            % Generate the kinetic part of the Hamiltonian for k(i_k,:)
            A0=spdiags(sum(abs(G+ones(nG,1)*k(i_k,:)).^2,2),0,nG,nG);
            % Calculate the eigenvalues
            E(i_k,:)=sort(real(eig(A0+A2)));
        end


        %%%% Plot the bandstructure in electronvolt
        plot((E-max(E(:,2)))*ryd_eV); axis([1 nk -6 6]); xlabel('k'); ylabel('E (eV)');
        %%%%
        
    elseif nargin==0 % Call the default values
        
        a0=5.43e-10/ryd_m;
        b=[ 0.53706 2.19104 2.05716 0.48716]; % Values from Zhang et al., Phys. Rev. B 48, 11204–11219 (1993) 
        
        n_atoms=2;
        tau=a0/8*[-1 1 ; -1 1 ; -1 1 ]; % Atom positions in the unit cell
        E_cut = 4.5;
        
        E = bulk_pseudo_william(a0, b, n_atoms, tau, E_cut);
        fprintf('The bandgap is %1.2feV\n',min(E(:,5)*ryd_eV)-max(E(:,2)*ryd_eV))
        
    end
return