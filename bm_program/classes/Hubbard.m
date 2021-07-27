classdef Hubbard < System_cl
    % creates general Hubbard Hamiltonian that conserves particle number and spin
    properties
        V_mat %
        U_mat %two particle interaction between different spins(site resolved)
        Hsp_mat %matrix containing single particle Hamiltonian in same spin / particle sector
        Geomet %Geometry describing the bonds, labelling the sites
        Part_hole_symmetric %describes if Hamiltonian is particle hole symmetric defined
        
    end
    
    
    methods
        function s = Hubbard(Geom, b_inp, U_inp, xi_inp, V_inp, Bas, part_hole_sym)
            %definition of parameters
            %Geom - defines Geometry of system (doesn't include any spin
            %   information but distance of atoms resp. arrangement of atoms
            %   who is neighbor, who not
            %b_inp - hopping, either a parameter for all, for each bond (as matrix / vector   )
            %   or even full matrix including on-site energies
            %U_inp - two particle interactions - where spin and particle
            %   number remain conserved for each site couples different spins
            %xi_inp - onsite energies, if not included already in b
            %V_inp - same as U but couples same spins
            %   comment: two particle operator for superconduction not
            %   implemented
            
            %check all input variables for different modes
            spin_symmetric = true;
            if nargin < 1 || isempty(Geom), Geom = Geometry(1); end
            if isinteger(Geom) && numel(Geom) == 1, Geom = Geometry(Geom); end %creates a linear chain with Geom sites
            if nargin < 2 || isempty(b_inp), b_inp = 0;
            elseif iscell(b_inp) && numel(b_inp) == 2, spin_symmetric = false; end %checks if hopping is set
            % if b is a cell - than all other parameters have to be set up as cells of two elements
            
            if nargin < 3 || isempty(U_inp), if spin_symmetric, U_inp = 0; else, U_inp =  {0,0}; end;
            elseif ~spin_symmetric && ~(iscell(U_inp) && numel(U_inp) == 1),error('U has to be set up as cell of one argument since it is not spin dependent defined! U_ij n_iup n_jdown'); end
            
            if nargin < 4 || isempty(xi_inp), if spin_symmetric, xi_inp = 0; else, xi_inp =  {0,0}; end;
            elseif ~spin_symmetric && ~(iscell(xi_inp)&& numel(xi_inp) == 2), error('xi is not set up for two spins');  end
            
            if nargin < 5 || isempty(V_inp), if spin_symmetric, V_inp = 0; else, V_inp = {0,0}; end;
            elseif ~spin_symmetric && ~(iscell(V_inp)&& numel(V_inp) == 2), error('V is not set up for two spins'); end
            
            if nargin < 6 || isempty(Bas), Bas = Basis(Geom.N_site); end
            
            if nargin < 7 || isempty(part_hole_sym), part_hole_sym = true; end
            
            %TODO: check through whole versions
            %check if Basis is complete
            if Bas.Numel_Bas ~= (2*Bas.N_spin)^Geom.N_site, warning('not using complete Basis set for Hubbard model'); end
            % SPIN: in general distinguish between different spins by
            % cells: b{1} (down) and b{2} (up)
            
            %initialize Table for Basis
            if Bas.N_spin == 2
                spin_used = 'S';
            else
                spin_used = [];
            end
            
            
            Tab = Table(Bas,['N', spin_used]);
            
            
            % cases for b: b is a number, b is a matrix of size
            % (N_site,N_site), b is a vector of Geom.bond_size, b is empty
            % (no hopping)
            Hsp = zeros(Geom.N_site);
            
            for k = 1:(~spin_symmetric+1)
                if ~spin_symmetric
                    
                    b = b_inp{k};
                    xi = xi_inp{k};
                    U = U_inp{1}; %U spin independent defined!
                    V = V_inp{k};
                    
                    %if ~isempty(U_inp{2}) && any(any(U_inp{1} ~= U_inp{2})), warning('U is not spin dependent defined, it is not symmetric!!! second cell discarded!'); end
                else
                    b = b_inp;
                    xi = xi_inp;
                    U = U_inp;
                    V = V_inp;
                end
                
                on_site_energy = true;
                if numel(b) == 1    % same hopping for all
                    Hsp = Geom.Bonds * b; %if Geom.bonds is not empty - but shouldn't be...
                    
                elseif numel(b) == Geom.N_bonds/2 %individual hopping for every bond ~= 0 (symmetric hopping)
                    %TODO: question about not symmetric hopping??? would violate hermiticity of Hamiltonian
                    b = b(:);
                    Hsp(Geom.Bond_index > 0)= Geom.Bonds(Geom.Bond_index > 0).*b(Geom.Bond_index(Geom.Bond_index>0));
                    Hsp(Geom.Bond_index < 0)= Geom.Bonds(Geom.Bond_index < 0).*b(-Geom.Bond_index(Geom.Bond_index<0));
                    %TODO: check if this works
                elseif all(size(b) == Geom.N_site) %
                    if all(diag(b) == 0)    % define full hopping (Geom.bond not needed / overruled) bonds can be set to 0
                        Hsp = b; %on site energy added later
                    else
                        %Matrix containing more entries than bonds - includes onsite energy
                        Hsp = b;
                        on_site_energy = false;
                    end
                else, error('size of b does not fit N_site')
                end
                
                if on_site_energy == false
                    warning('on-site energy already defined in b')
                    if spin_symmetric
                        xi_inp = [];
                    else
                        xi_inp = {[],[]};
                    end
                end
                if on_site_energy  %case when an explicit on-site energy is used
                    if numel(xi) == 1
                        Hsp(eye(Geom.N_site) == 1) = xi;
                    elseif numel(xi) == Geom.N_site
                        Hsp(eye(Geom.N_site) == 1) = xi(:);
                    elseif all(size(xi) == size(Hsp))
                        Hsp(eye(Geom.N_site) == 1) = diag(xi);
                    else, error('size of xi does not fit N_site')
                    end
                end
                
                %necessary for input already??? -> at least warning
                Hsp_test = (Hsp +Hsp')/2;
                if any(Hsp(:) ~= Hsp_test(:))
                    warning('input parameters b and xi not apt to make H self-adjoint')
                    %Hsp = Hsp_test;
                    % update input, since input is corrected here to make Hsp self-adjoint
                    if spin_symmetric
                        b_inp = Hsp;
                    else
                        b_inp{k} = Hsp;
                    end
                end
                
                
                
                
                
                
                %U - note: next nearest neighbour influence of two particle
                %interaction not implemented via Geom.bonds since no extra
                %parameter for this - best way is to create matrix for U
                %and V
                
                if numel(U) == 1 %homogeneous U
                    Umat = eye(Geom.N_site)*U;
                elseif numel(U) == Geom.N_site %U different for each site
                    Umat = diag(U);
                elseif size(U,1) == Geom.N_site && size(U,2) == Geom.N_site  % more electron-electron interaction of different spin
                    Umat = U;
                else, error('size of U does not fit N_site')
                end
                if Bas.N_spin == 1
                    if nnz(Umat) ~= numel(Umat), warning('for spinless fermions U set to zero'); end
                    Umat = Umat * 0;
                end
                %usage of U - U does not have a spin index, since it describes the interaction between opposite
                %spins - the definition is the following: U_{ij} n_{i\uparrow} n_{j\downarrow}
                % U is not symmetric since U_{12} is different than U_{21}:   n_{1up} n_{2down} vs n_{2up} n_{1down}
                
                Vmat = zeros(Geom.N_site,Geom.N_site);
                if numel(V) == 1 %homogeneous V for all bonds
                    Vmat(Geom.Bonds ~= 0) = V;
                elseif numel(V) == Geom.N_bonds/2 %V different for each bond
                    V = V(:);
                    Vmat(Geom.Bond_index > 0)= V(Geom.Bond_index(Geom.Bond_index>0)); % Geom.Bonds(Geom.Bond_index > 0) -> TODO: maybe introduce an own Geometry to describe electron electron interactions
                    Vmat(Geom.Bond_index < 0)= V(-Geom.Bond_index(Geom.Bond_index<0));
                elseif size(V,1) == Geom.N_site && size(V,2) == Geom.N_site  % more electron-electron interaction of different spin
                    if any(diag(V) ~= 0)
                        error('violates Pauli principle - there is no interation of the same spin at the same site')
                    end
                    Vmat = V;
                else, error('size of V does not fit N_site')
                end
                if any(any(Vmat ~= Vmat.'))
                    warning('interaction V is bidirectional - there should not be a difference!')
                end
                % since V is bidirectional lower triangular part is deleted
                % n_{i up} n_{j up} = n{j up} n_{i up}
                Vmat(tril(true(size(Vmat)),-1)) = 0;
                
                
                if ~spin_symmetric
                    st.U_mat{1} = Umat; %since U is not spin dependent!!!
                    st.Hsp_mat{k} = Hsp;
                    st.V_mat{k} = Vmat;
                else
                    st.U_mat = Umat;
                    st.Hsp_mat = Hsp;
                    st.V_mat = Vmat;
                end
            end
            
            %==============================================================
            %--------------------------------------------------------------
            % create many body Hamiltonian
            %--------------------------------------------------------------
            
            % TODO: if is particle/spin seperated, create a procedure that works just in each sector
            % especially work in spin / particle sector for two particle situations
            
            
            
            % first: the question is here whether to work in the
            % individual particle/spin sectors, since transformation of
            % all basis vectors to binary could become cost
            % intensive,
            
            % second work with unique spin up binary part and unique spin
            % down binary part : kron(Ham,Ham)
            %sectors - both already in sectors:
            %spin up
            bas_bin = Bas.f_bin();
            
            if lower(Bas.Smaller_Binary_Sector_Used_For(1)) == 'd'
                % indicates which spin sector is on the right side 0101|1001
                %spin_up_factor gives digits to shift from one spin sector
                %to the other in binary basis (from the left)
                spin_up_shift = 0;
            else, spin_up_shift = Bas.N_site;
            end
            spin_down_shift = Bas.N_site - spin_up_shift;
            
            
            % TODO: verify this rough estimation of 1 percent for number of nonzero elements
            %Hamiltonian_mb = spalloc(Bas.Numel_Bas,Bas.Numel_Bas, );
            [I1, I2, W] = deal(zeros(Bas.Numel_Bas*Bas.N_spin*Bas.N_site,1));
            sparse_index = 0;
            
            % ----------Hsp part (onsite energy and hopping)----------------------
            %spin up
            % if spin symmetric all is done for both spins, so \sum_{\sigma, ij} t_ij
            for spin = 1:Bas.N_spin
                % 1...spin down, 2...spin up
                if ~spin_symmetric
                    Hsp =  st.Hsp_mat{spin};
                    Umat = st.U_mat{1};  %since U is not spin dependent!
                    Vmat = st.V_mat{spin};
                else
                    Hsp = st.Hsp_mat;
                    Umat = st.U_mat;
                    Vmat = st.V_mat;
                end
                if spin == 1, spin_shift = spin_down_shift;
                else, spin_shift = spin_up_shift;
                end
                
                if Bas.N_spin == 1 % spinless fermions
                    spin_shift = 0;
                end
                
                
                for ind_i = spin_shift + (1:Geom.N_site)
                    for ind_j = spin_shift + (1:Geom.N_site)
                        if ~isdeployed()
                            disp(['i: ', num2str(ind_i), ' to j: ', num2str(ind_j)])
                        end
                        
                        if Hsp(ind_i-spin_shift, ind_j-spin_shift) ~= 0
                            if ind_i == ind_j
                                %L2 = bas_bin(:,ind_i) == 1;
                                %LL = diag(L2) == 1;
                                
                                L2 = find(bas_bin(:,ind_i) == 1);
                                % get diagonal index
                                %LL = (L2-1)*(Bas.Numel_Bas+1) +1;
                                %Hamiltonian_mb(LL) = Hamiltonian_mb(LL) + Hsp(ind_i-spin_shift, ind_i-spin_shift);
                                
                                % on site energy
                                
                                
                               [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, L2, L2, Hsp(ind_i-spin_shift, ind_i-spin_shift) * ones(size(L2)));
                              
                                
%                                 sparse_index = max(sparse_index) + (1:numel(L2));
%                                 I1(sparse_index) = L2;
%                                 I2(sparse_index) = L2;
%                                 W(sparse_index) =  Hsp(ind_i-spin_shift, ind_i-spin_shift) * ones(size(L2));
%                                 
                            else
                                
                                L2 = find((bas_bin(:,ind_i) == 0) & bas_bin(:, ind_j) == 1); %ket vectors, columns in H
                                if isempty(L2)
                                    warning('used basis set not complete - so some hoppings not possible')
                                else
                                    num_ket = numel(L2);
                                    % for Fermi sign extra minus if j < i (has to pass), since we always have: c_i^dagger c_j
                                    F2 = (-1).^(sum([zeros(num_ket,1), bas_bin(L2, 1:(ind_i-1))],2) + sum([zeros(num_ket,1), bas_bin(L2,1:ind_j-1)],2) + (ind_j < ind_i));
                                    Bas_new = Bas.Bas_Dez(L2)+2^(Bas.N_spin*Geom.N_site-ind_i) - 2^(Bas.N_spin*Geom.N_site - ind_j);
                                    [~,L1] = ismember(Bas_new,Bas.Bas_Dez); %optimization for search function if part/spin sectors are preserverd!!!
                                    if all(~L1)
                                        %warning('break - used basis set not complete - so some hoppings not possible')
                                    else
                                        if any(~L1)
                                            %warning('used basis set not complete - so some hoppings not possible')
                                            L1 = L1(L1~=0);
                                            L2 = L2(L1~=0);
                                            F2 = F2(L1~=0);
                                        end
                                        
                                        [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, L1, L2, Hsp(ind_i-spin_shift, ind_j-spin_shift) .*F2);
%                                         sparse_index = max(sparse_index) + (1:numel(L1));
%                                         I1(sparse_index) =  L1;
%                                         I2(sparse_index) =  L2;
%                                         W(sparse_index) =  Hsp(ind_i-spin_shift, ind_j-spin_shift) .*F2;
                                        
                                        
                                        %LL = sub2ind([Bas.Numel_Bas,Bas.Numel_Bas], L1,L2); %TODO: test if sub2ind is slow or not (could be, because of argument check)
                                        %hopping
                                        %Hamiltonian_mb(LL) = Hamiltonian_mb(LL) + Hsp(ind_i-spin_shift, ind_j-spin_shift) .*F2;
                                    end
                                end
                                
                            end
                        end
                        
                        %Umat
                        if Umat(ind_i-spin_shift, ind_j-spin_shift) ~= 0
                            %U just runs through for one spin
                            
                            if spin == 1 %indices are spin_down, ind_i must be corrected to up
                                % spin == 1
                                % ind_i ... spin_up
                                % ind_j ... spin_down
                                % logical vectors
                                Lleft = bas_bin(:,ind_i-spin_down_shift + spin_up_shift) == 1;
                                Lright = bas_bin(:,ind_j) == 1;
                                LU = find(Lleft & Lright);
                                if isempty(LU)
                                    warning('used basis set not complete - so some U interactions not possible')
                                end
                                
                                % get the diagonal indices
                                %LL = (find(LU)-1)*(Bas.Numel_Bas+1) +1;
                                %LLu = (find(Lleft)-1)*(Bas.Numel_Bas+1) +1;
                                %LLd = (find(Lright)-1)*(Bas.Numel_Bas+1) +1;
                                %Hamiltonian_mb(LL) = Hamiltonian_mb(LL) + Umat(ind_i-spin_shift, ind_j-spin_shift);
                                
                                
                                LR = find(Lright);
                                LL = find(Lleft);
                                
                                [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, LU, LU, Umat(ind_i-spin_shift, ind_j-spin_shift) * ones(size(LU)));

%                                 sparse_index = max(sparse_index) + (1:numel(LU));
%                                 I1(sparse_index) = LU;
%                                 I2(sparse_index) = LU;
%                                 W(sparse_index)  = Umat(ind_i-spin_shift, ind_j-spin_shift) * ones(size(LU));
%                                 
                                
                                
                                % make system particle hole symmetric (maybe
                                if part_hole_sym
                                    %Hamiltonian_mb(LLu) = Hamiltonian_mb(LLu) - 1/2*Umat(ind_i-spin_shift, ind_j-spin_shift);
                                    %Hamiltonian_mb(LLd) = Hamiltonian_mb(LLd) - 1/2*Umat(ind_i-spin_shift, ind_j-spin_shift);
                                    [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, LL, LL, - 1/2*Umat(ind_i-spin_shift, ind_j-spin_shift) * ones(size(LL)));
                                    [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, LR, LR, - 1/2*Umat(ind_i-spin_shift, ind_j-spin_shift) * ones(size(LR)));
                                    
%                                     sparse_index = max(sparse_index) + (1:(numel(LL) + numel(LR)));
%                                     I1(sparse_index) = [LL; LR];
%                                     I2(sparse_index) = [LL; LR];
%                                     W(sparse_index)  = [- 1/2*Umat(ind_i-spin_shift, ind_j-spin_shift) * ones(size(LL)); - 1/2*Umat(ind_i-spin_shift, ind_j-spin_shift) * ones(size(LR))];
                                end
                                
                                %spin == 2 not used, since it is a repitition
                            end
                            
                        end
                        
                        %Vmat!!!!!!!!
                        % Vmat is nonzero just for upper triangle since n_is n_js = n_js n_is
                        % though performed for both spins. diag(Vmat) == 0
                        if Vmat(ind_i-spin_shift, ind_j - spin_shift) ~= 0
                            L_first = bas_bin(:,ind_i) == 1;
                            L_second = bas_bin(:,ind_j) == 1;
                            LV = find(L_first & L_second);
                            if isempty(LV)
                                    warning('used basis set not complete - so some V interactions not possible')
                            end
                            %LL = (find(LV)-1)*(Bas.Numel_Bas+1) +1;
                            %LLfirst = (find(L_first)-1)*(Bas.Numel_Bas+1) +1;
                            %LLsecond = (find(L_second)-1)*(Bas.Numel_Bas+1) +1;
                            %Hamiltonian_mb(LL) = Hamiltonian_mb(LL) + Vmat(ind_i-spin_shift, ind_j - spin_shift);
                            
                            
                            
                            [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, LV, LV, Vmat(ind_i-spin_shift, ind_j - spin_shift) * ones(size(LV)));
%                             sparse_index = max(sparse_index) + (1:numel(LV));
%                             I1(sparse_index) =  LV;
%                             I2(sparse_index) =  LV;
%                             W(sparse_index)  = Vmat(ind_i-spin_shift, ind_j - spin_shift) * ones(size(LV));
                            
                            
                            
                            % TODO: check if not good idea for U acting on neighboring sites
                            % make system particle hole symmetric
                            if part_hole_sym
                                %Hamiltonian_mb(LLfirst) = Hamiltonian_mb(LLfirst) - 1/2*Vmat(ind_i-spin_shift, ind_j - spin_shift);
                                %Hamiltonian_mb(LLsecond) = Hamiltonian_mb(LLsecond) - 1/2*Vmat(ind_i-spin_shift, ind_j - spin_shift);
                                L1 = find(L_first);
                                L2 = find(L_second);
                                
                                [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, L1, L1, - 1/2*Vmat(ind_i-spin_shift, ind_j - spin_shift) * ones(size(L1)));
                                [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, L2, L2, - 1/2*Vmat(ind_i-spin_shift, ind_j - spin_shift) * ones(size(L2)));
                                
%                                 sparse_index = max(sparse_index) + (1:(numel(L1) + numel(L2)));
%                                 I1(sparse_index)  = [L1; L2];
%                                 I2(sparse_index)  = [L1; L2];
%                                 W(sparse_index) = [- 1/2*Vmat(ind_i-spin_shift, ind_j - spin_shift) * ones(size(L1)); - 1/2*Vmat(ind_i-spin_shift, ind_j - spin_shift) * ones(size(L2))];
%                                 
                            end
                            
                        end
                        
                    end
                end
                
                
                
                
                %-----------two particle Hamiltonian---------------
                % TODO:
                
            end
            
            % check double entries not necessary
            if ~(numel(I1) > max(sparse_index) )
                warning(['estimation of nonzero elements (', num2str(Bas.Numel_Bas*Bas.N_spin*Bas.N_site), ...
                    ') was too weak, nnz= ', num2str(max(sparse_index))]);
            else
                I1 = I1(1:max(sparse_index));
                I2 = I2(1:max(sparse_index));
                W = W(1:max(sparse_index));
            end
            
            
            %addition of same indices...
            
%             A = spalloc(Bas.Numel_Bas, Bas.Numel_Bas,Bas.Numel_Bas*30 )
%             A = sparse(I1,I2,W,Bas.Numel_Bas, Bas.Numel_Bas);
%             [~,ix] = sort(I1);
%             ind2 = uint64(sub2ind([Bas.Numel_Bas, Bas.Numel_Bas], I1, I2));
%             
%             ind2(1:30)
%             
            
            
%             addpath GPU
%             comment = 'three band hubbard model CuO, N_site = 9, N_e = 45/54 electrons, S = 9';
%             mmwrite_sparse(['/temp/dorn/run/CuO/hubbard3x3_S9.mtx'], I1,I2,W,[Bas.Numel_Bas, Bas.Numel_Bas],comment,'hermitian')
%             

            
            Hamiltonian_mb = sparse(I1,I2,W,Bas.Numel_Bas,Bas.Numel_Bas);
            
            
            
            %TODO:  pcolor(full(Hamiltonian_mb))
            
            % create
            
            
            
            
            
            
            % if spin_symmetric %assign second block
            Model_name = 'hubbard';
            
            model_name = {Model_name};
            geom_struc = {Geom.description };
            if iscell(b_inp), hopping = b_inp; else, hopping = {b_inp}; end
            if iscell(U_inp), interaction = U_inp; else, interaction = {U_inp}; end
            if iscell(xi_inp), on_site_energy = xi_inp; else, on_site_energy = {xi_inp}; end
            if iscell(V_inp), same_spin_inter = V_inp; else, same_spin_inter = {V_inp}; end
            N_site = Geom.N_site;
            N_spin = Bas.N_spin;
            particle_hole_symmetric = part_hole_sym;
            
            % Call asset constructor
            
            parameter_list = table(model_name, geom_struc,N_site,spin_symmetric,  hopping,on_site_energy, interaction,  same_spin_inter, N_spin, particle_hole_symmetric);
            
            s@System_cl(parameter_list, Model_name, Bas, Hamiltonian_mb, Tab);
            s.Hsp_mat = st.Hsp_mat;
            s.U_mat = st.U_mat;
            s.V_mat = st.V_mat;
            s.Geomet = Geom;
            s.Part_hole_symmetric = part_hole_sym;
        end
        
        function s_print(obj)
            s_print@System_cl(obj);
            
            
            
            
        end
        
        function [hamiltonian, parameters] = s_hamiltonian(obj)
            %hamiltonian = ['$$\sum_{i=1}^', num2str(obj.Parameter_list.N_site), 'c_i^\dag c_j t_{ij} + U \sum_i n_i n_{i+1}$$'];
            
            if obj.Parameter_list.N_spin == 1 %spin symbol always used in sums if there is a spin
                spin_symbol = '';
                spin_sum= '';
            else
                spin_symbol = '\sigma';
                spin_sum = '\sigma, ';
            end
            
            if obj.Parameter_list.spin_symmetric == false %param_spin always used if parameters are spin_dependent
                param_spin = '\sigma';
            else
                param_spin = '';
            end
            
            on_site_param = obj.Parameter_list.on_site_energy{1};
            
            if numel(on_site_param) == 1
                on_site = ['\xi_{',param_spin, '}'];
            else
                on_site = ['\xi_{i', param_spin, '}'];
            end
            if isempty(on_site_param)
                hamiltonian_on_site = '';
            else
                hamiltonian_on_site = ['\sum_{', spin_sum, 'i=1}^N', on_site, ' n_{i', spin_symbol, '}+'];
            end
            
            
            
            
            
            hopping_param = obj.Parameter_list.hopping{1};
            
            hopping_sum_start = [spin_sum, 'i = 1'];
            hopping_sum_boundary = '';
            sec_hopp_op_ind = 'j';
            if numel(hopping_param) == 1
                hopping = ['t_{', param_spin, '}'];
            elseif numel(hopping_param) == obj.Geomet.N_bonds/2
                hopping = ['t_{i', param_spin, '}'];
                if strcmp(obj.Parameter_list.geom_struc,'periodic cycle')
                    hopping_sum_boundary = '^N';
                    sec_hopp_op_ind = 'i+1';
                elseif strcmp(obj.Parameter_list.geom_struc, 'linear chain')
                    hopping_sum_boundary = '^{N-1}';
                    sec_hopp_op_ind = 'i+1';
                end
            else
                hopping = ['t_{ij', param_spin, '}'];
                hopping_sum_start = [spin_sum, 'ij'];
                hopping_sum_boundary = '^N';
            end
            
            
            
            
            hamiltonian_hopping = ['\sum_{', hopping_sum_start, '}', hopping_sum_boundary, hopping '(c_{i', spin_symbol, '}^\dag ', 'c_{', sec_hopp_op_ind, spin_symbol, '} + h.c.)'];
            % ij, <ij> i (i+1)  c_i^\dag c_j  c_i^\dag c_{i+1}
            
            %interaction
            
            if obj.Parameter_list.N_spin > 1
                interaction_param = obj.Parameter_list.interaction{1};
                interaction_sum_boundary = '^N';
                sec_interaction_op_ind = 'i'; %no summation over spin, order {up,i} {down,j} fixed
                interaction_sum_start = 'i = 1';
                if numel(interaction_param) == 1
                    interaction = 'U';
                    
                elseif numel(interaction_param) == obj.Parameter_list.N_site
                    interaction = 'U_{i}';
                else
                    interaction = 'U_{i j}';
                    interaction_sum_start = 'ij'; %no summation over spin, order up, down fixed
                    sec_interaction_op_ind = 'j';
                    interaction_sum_boundary = '^N';
                end
                
                part_hole_begin_interaction = '';
                part_hole_end_interaction = '';
                if obj.Part_hole_symmetric
                    part_hole_begin_interaction = '\left(';
                    part_hole_end_interaction = [' - \frac{n_{i\uparrow} + n_{', sec_interaction_op_ind, '\downarrow}}{2}\right)'];
                end
                
                hamiltonian_interaction = ['+ \sum_{' interaction_sum_start, '}', interaction_sum_boundary, interaction, part_hole_begin_interaction, 'n_{i\uparrow} n_{', sec_interaction_op_ind, '\downarrow}', part_hole_end_interaction];
                %i, ij, i\sigma, ij\sigma   (U_i,U) n_{i\uparrow} n_{i\downarrow}, U_{ij} n_{i\uparrow} n_{j \downarrow}, (U_{i\simga}, U_\sigma) n_{i\sigma} n_{i\overline{\sigma}, U_{ij\sigma} n_{i\sigma} n_{j\overline{\sigma}
            else
                hamiltonian_interaction = '';
            end
            % + particle hole symmetric!!!
            % TODO: write a list with all cases and look how to map which one to which
            %same_spin_inter_param
            same_spin_inter_param = obj.Parameter_list.same_spin_inter{1};
            ssi_sum_start = [spin_sum, 'i = 1'];
            sec_ssi_op_ind = 'j';
            ssi_sum_boundary = '';
            if numel(same_spin_inter_param) == 1
                ssi = ['V_{', param_spin, '}'];
            elseif numel(same_spin_inter_param) == obj.Geomet.N_bonds/2
                ssi = ['V_{i', param_spin, '}'];
                if strcmp(obj.Parameter_list.geom_struc,'periodic cycle')
                    ssi_sum_boundary = '^N';
                    sec_ssi_op_ind = 'i+1';
                elseif strcmp(obj.Parameter_list.geom_struc, 'linear chain')
                    ssi_sum_boundary = '^{N-1}';
                    sec_ssi_op_ind = 'i+1';
                end
            else
                ssi = ['V_{ij', param_spin, '}'];
                ssi_sum_start = [spin_sum, 'i > j'];
                ssi_sum_boundary = '^N';
            end
            
            part_hole_begin_ssi = '';
            part_hole_end_ssi = '';
            if obj.Part_hole_symmetric
                part_hole_begin_ssi = '\left(';
                part_hole_end_ssi = [' - \frac{n_{i', spin_symbol, '} + n_{', sec_ssi_op_ind, spin_symbol,'}}{2}\right)'];
            end
            
            
            hamiltonian_same_spin_inter = ['+ \sum_{' ssi_sum_start, '}', ssi_sum_boundary, ssi, part_hole_begin_ssi, 'n_{i',spin_symbol, '} n_{', sec_ssi_op_ind, spin_symbol, '}', part_hole_end_ssi];
            
            
            hamiltonian = ['$$ H_S = ', hamiltonian_on_site, hamiltonian_hopping, hamiltonian_interaction, hamiltonian_same_spin_inter, ' $$'];
            
            if obj.Parameter_list.spin_symmetric
                k_max = 1;
            else
                k_max = 2;
            end
            
            for k = 1:k_max
                
                hopping_param = obj.Parameter_list.hopping{k};
                interaction_param = obj.Parameter_list.interaction{1}; %since U is never spin dependent
                same_spin_inter_param = obj.Parameter_list.same_spin_inter{k};
                on_site_param = obj.Parameter_list.on_site_energy{k};
                
                if all(size(hopping_param) > 1)
                    hopping_parameter_text = 'matrix';
                else
                    hopping_parameter_text = num2str(transpose(hopping_param(:)));
                end
                
                if all(size(interaction_param) > 1)
                    interaction_parameter_text = 'matrix';
                else
                    interaction_parameter_text = num2str(transpose(interaction_param(:)));
                end
                
                if all(size(same_spin_inter_param) == obj.Parameter_list.N_site)
                    ssi_parameter_text = 'matrix';
                else
                    ssi_parameter_text = num2str(transpose(same_spin_inter_param(:)));
                end
                
                on_site_parameter_text = num2str(transpose(on_site_param(:)));
                
                parameters(1:2,k) = {['$$ ', on_site, '$$: ', on_site_parameter_text]; ...
                    ['$$ ', hopping, '$$: ', hopping_parameter_text, ]};
                if obj.Parameter_list.N_spin == 1
                    parameters(3,k) = { ['$$ ', ssi, '$$: ', ssi_parameter_text]};
                else
                    parameters(3:4,k) = {  ['$$ ', interaction, '$$: ', interaction_parameter_text]; ...
                        ['$$ ', ssi, '$$: ', ssi_parameter_text]};
                end
                
                
            end
            if k_max == 2
                parameters = [{'\underline{spin down}', '\underline{spin up}'}; parameters];
            end
        end
        
    end
end