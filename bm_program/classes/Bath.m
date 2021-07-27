classdef Bath
    %author: Gerhard Dorn
    %date: October 19
    %use: describes the non-interacting bath including orbitals used in the master equation approach, bath is described by Temperature (beta), Energy
    %and bias dependence, filling type (dagger), Gamma / hopping inside is rather fixed / constant
    % 
    %features: 
    %          handle:
    %             - uses precalculated values to fasten up the access to them (uses interp2 / interp3 - if
    %             temperature dependent
    %             - expression for Greens function
    %             - saves bath type, and typical density of states used
    %           
    %             
    %          visualize:
    %             - print / visualize bath / energies in relation to system
    %          initialize:
    %             - performs all integrals to calculate Bath_corr_func
    %             
    
    % TODO: check implications of a more general bath (at the moment we use particle number momentum and spin
    % conservation (not true for superconductors) and different baths are not correlated
      
    properties
        
        N_spin          %NEW
        Numel_Baths 
        Beta
        Mu
        Spin_conserved  %NEW indicates whether spin is conserved or not (spinless, no mixing of spin up and down)
        Spin_symmetric  %NEW indicates whether spin is symmetric (just one sector updown stored, spinless)
        Numel_Orbitals
        
        Omega      % always defined in the beginning
        GF  % Green functions describing the bath
        
        Fermi =  @(x,mu,beta,dag) 1./(exp(dag*beta * (x - mu)) + 1);
        Shift_Baths     % True of false if shift of the chemical potential to the Green's function is performed or not.
                
        Integrator      % defines which algorithm is used for the integration
        Bath_Corr       % used to create a table to store already calculated correlation functions
        Bath_Calc       % used to create a table to store already calculated correlation functions
        En_Range        % used to create a table to store already calculated correlation functions
        Mu_Range        % used to create a table to store already calculated correlation functions
        Beta_Range      % used to create a table to store already calculated correlation functions
        
        

        
        
        
        
    end
    
    methods
        
        function obj = Bath(GR, grid, shift_baths, integrator, N_spin, Numel_orbitals, Spin_conserved, Spin_symmetric, Beta, Mu)
            if nargin < 9 || isempty(Beta), Beta = 20*ones(size(GR,1),1); end
            if nargin < 10 || isempty(Mu), Mu = zeros(size(GR,1),1); end
            
            
            obj.N_spin = N_spin;      
            obj.Numel_Baths = size(GR,1);
            obj.Beta = Beta;
            obj.Mu = Mu;
            obj.Spin_conserved = Spin_conserved; %NEW indicates whether spin is conserved or not (spinless, no mixing of spin up and down)
            obj.Spin_symmetric = Spin_symmetric; %NEW indicates whether spin is symmetric (just one sector updown stored, spinless)
            obj.Numel_Orbitals = Numel_orbitals;

            
            if ~isa(grid, 'Grid')
                obj.Omega = Grid(grid);
            else
                obj.Omega = grid;
            end
            
            obj.GF = GR;
            

            if nargin < 3 || isempty(shift_baths), shift_baths = false(obj.Numel_Baths,1);end
            obj.Shift_Baths = shift_baths;
            if nargin < 4 || isempty(integrator), integrator = 'Fast'; end
            obj.Integrator = integrator;
            
            
            %TODO: treat also noncell GR input!!!!
            
               %NEW
            
            
%             if obj.N_spin == 1
%                 obj.Spin_fun = @(x) 1; %TODO: rethink Spin_fun!!! - when and for what is it used? addressing the right sector in GF?
%             elseif obj.N_spin == 2
%                 if obj.Spin_conserved && obj.Spin_symmetric
%                     obj.Spin_fun = @(x) 1;
%                 else
%                     obj.Spin_fun = @(x) x/2+1.5;
%                 end
%             else
%                 error('unusual number of spins! Not implemented!')
%             end
            
         
            %Two strategies to store the bath correlation values:
            % 1) as Category
            
            % 2) as sparse
            
            % question about preallocation, how many Delta E are needed?
            
            
            % Category used in a bit different way...
            % add single entries and delete single entries according to all Categories
            % ----> Slow 
            %TODO: think of a faster strategy (preallocated sparse matrix)
        
        
        end
        
        function Numel_Spins(obj, nargin)
            error('Numel_Spins in Bath class: to be deleted!')
        end
        
        function Spin_fun(obj, nargin)
            error('Spin_fun in Bath class: to be deleted!')
        end
        %
        function GF = f_GF_shifted(obj, bath_id, mu)
            %Define what this function is doing.
            %NOTE: used in CPT_current!
            
            if  numel(obj.GF{1}) < 2
                error('GF not defined (constant bath), CPT will not work.')
            end
                
            if strcmp(obj.Integrator, 'Constant') 
                warning('Do not use constant baths with CPT')
            end
            
            if ~obj.Shift_Baths(bath_id) || mu == 0 
                GF = obj.GF{bath_id, 1};
            else
                [nx,ny,nz] = size(obj.GF{bath_id, 1});
                GF_extended = zeros(nx,ny,nz+1);
                if mu > 0
                    grid_extended = [obj.Omega.Limits(1)- mu; obj.Omega.Points];
                    GF_extended(1:nx,1:ny,2:end) = obj.GF{bath_id, 1};
                    %extrapolate boundary point!
                elseif mu < 0
                    grid_extended = [obj.Omega.Points; obj.Omega.Limits(2) - mu]; %mu is negative
                    GF_extended(1:nx,1:ny,1:end-1) = obj.GF{bath_id, 1};
                    %extrapolate boundary point!
                end
                
                GF = reshape(spline(grid_extended, GF_extended, obj.Omega.Points- mu),nx,ny,nz);
            end
        end
            
        function C = f_gamma_mat(obj, omega, alpha, dag_mat)
            % calculates i(GR - GA) f^{\overline{s}}(s\omega)
            % four things to consider:
            % 1) Argument used is dag_mat .* omega so do not do it before!
            % 2) s of fermi function is given explicitly as dag_mat (f^{\overline{s}})
            % 3) shift is done as with f_gamma take original GF and shift the corresponding data points by -
            % shift\mu
            % 4) permutation is included with dag_mat
            % works for orbitas, generalization to non-equilibrium bath not considered yet.
            
            
            GFAK_s_plus = 1i*(obj.GF{alpha,1}(:,:,:) - conj(permute(obj.GF{alpha,1}, [2,1,3]))) .* ...
                repmat(permute(1-obj.Fermi(obj.Omega.Points, obj.Mu(alpha)+obj.Shift_Baths(alpha)*obj.Mu(alpha), obj.Beta(alpha), 1), [2,3,1]), size(obj.GF{alpha,1},1), size(obj.GF{alpha,1},2));
            GFAK_s_minus = conj(1i*(obj.GF{alpha,1}(:,:,:) - conj(permute(obj.GF{alpha,1}, [2,1,3])))) .* ...
                repmat(permute(obj.Fermi(obj.Omega.Points, obj.Mu(alpha)+obj.Shift_Baths(alpha)*obj.Mu(alpha), obj.Beta(alpha), 1), [2,3,1]), size(obj.GF{alpha,1},1), size(obj.GF{alpha,1},2));
            
            L_s_plus = dag_mat == 1;
            L_s_minus = dag_mat == -1;
            C = zeros( [size(omega,1).* size(omega,2), obj.Numel_Orbitals(alpha), obj.Numel_Orbitals(alpha)]);
            
            C(L_s_plus(:),:,:) = interp1(obj.Omega.Points  -obj.Shift_Baths(alpha)*obj.Mu(alpha), permute(GFAK_s_plus, [3,1,2]), dag_mat(L_s_plus).*omega(L_s_plus), 'spline');
            C(L_s_minus(:),:,:) = interp1(obj.Omega.Points  -obj.Shift_Baths(alpha)*obj.Mu(alpha), permute(GFAK_s_minus, [3,1,2]), dag_mat(L_s_minus).*omega(L_s_minus), 'spline');
            C = reshape(C, [size(omega,1), size(omega,2), obj.Numel_Orbitals(alpha), obj.Numel_Orbitals(alpha)]);
            
        end
            
            
            
        function [C, Bat, precalc_warning] = f_gamma(obj, omega, mu_alpha, beta_alpha, alpha, c_spin, dag)
            
            % gamma = i/2 sum    GR-GA + sGK (dag\omega)   |_ \pi(dag)(lk)    ...   exchange l and k if dag is negative
            % in equilibrium K is given by (GR-GA) (1-2fermi)
            
            % sum = \sum_\alpha, kl   V_\alpha,k\mu^s V_\alpha,l\kappa^\overline(s)
            
            % in bath equilibrium!
            %NOTE: Bath Shift: Shift G^R and G^A but not  Fermi function! (G_K)
            %G_old(x) -> G_new(x) := G_old(x-mu) (if old peak at zero, new peak at mu), recorrect Fermi function in
            %GFAK: + mu, in interpolation substract: x-mu
            
            GFAK = 1i/2*(obj.GF{alpha,1}(:,:,:) - conj(permute(obj.GF{alpha,1}, [2,1,3]))) .* ...
                repmat(permute((1+dag*(1-2*obj.Fermi(obj.Omega.Points, mu_alpha+obj.Shift_Baths(alpha)*mu_alpha, beta_alpha, 1))), [2,3,1]), size(obj.GF{alpha,1},1), size(obj.GF{alpha,1},2));
            C = permute(interp1(obj.Omega.Points  -obj.Shift_Baths(alpha)*mu_alpha, permute(GFAK, [3,1,2]), dag*omega, 'spline'), [2,3,1]);
            C(isnan(C)) = 0;
            %transpose C if dag is negative (\pi(dag) (lk)
            
            if dag == -1
                C = permute(C, [2,1,3]);
            end
            
            %the question is, if we should use precalculation? - maybe not...
            Bat = obj;
            precalc_warning = false;
%             
%             pp = spline(obj.Omega.Points, squeeze(obj.GF{alpha, obj.Spin_fun(c_spin)}(row, column, :)) - conj(squeeze(obj.GF{alpha, obj.Spin_fun(c_spin)}(column, row, :))));
%             
%             for en_ind = 1:numel(energ_inp)
%                 energy = energ_inp(en_ind);
%                 
%                 func = @(x) - ppval(pp,x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(x,mu_alpha, beta_alpha, dag)./(2*pi *(dag * x + energy + 1i*obj.Omega.Imag_Infinitesimal));
%                 
%                 %C(-E)
%                 % NOTE: Bath shift by mu_alpha optionally included
%                 C(row, column, en_ind) = integral(@(x) - ppval(pp,x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(x,mu_alpha, beta_alpha, dag)./(2*pi *(dag * x + energy + 1i*obj.Omega.Imag_Infinitesimal)),...
%                     obj.Omega.Limits(1), obj.Omega.Limits(2));
%     %              ZZ(en_ind) = integral(@(x) - ppval(pp,x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(x,mu_alpha, beta_alpha, dag)./(2*pi *(dag * x + energy)),...
%                             obj.Omega.Limits(1), -dag * energy - epsi) + ...
%                             integral(@(x) - ppval(pp,x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(x,mu_alpha, beta_alpha, dag)./(2*pi *(dag * x + energy)),...
%                             -dag * energy + epsi, obj.Omega.Limits(2)) + ...
%                             integral(@(x) ((- ppval(pp,-dag*energy + x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(-dag*energy + x,mu_alpha, beta_alpha, dag))   - ...
%                                           (- ppval(pp,-dag*energy - x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(-dag*energy - x,mu_alpha, beta_alpha, dag))) ./(x*dag), 0, epsi);
                        
            
        end
        
        function [C,Bat] = f_sigma_mat(obj, energies, alpha, dag_mat)
            
            % to get W multiply with -i/2!!!
            % same as f_sigma!
            
            %shift baths works the following way:
            % 1) either: [omega, GF(omega) * fermi(omega, mu)] -> [omega, GF(omega+shift) * fermi(omega,mu)]
            % 2) or:     [omega, GF(omega) * fermi(omega, mu)] -> [omega - shift, GF(omega) * fermi(omega - shift,mu)]
            
            % here: 2
            l = 1; p= 1;
            func_p =   1/(pi)*squeeze(obj.GF{alpha}(l,p,:) - conj(obj.GF{alpha}(p,l,:)) ).*obj.Fermi(obj.Omega.Points-obj.Shift_Baths(alpha)*obj.Mu(alpha), obj.Mu(alpha),obj.Beta(alpha),  -1);
            func_m =   1/(pi)*conj(squeeze(obj.GF{alpha}(l,p,:)- conj(obj.GF{alpha}(p,l,:)))).*obj.Fermi(obj.Omega.Points-obj.Shift_Baths(alpha)*obj.Mu(alpha), obj.Mu(alpha), obj.Beta(alpha), 1);
            
            %interpolate to equidistant grid!
            % find support:
            left_index = min(find(cumsum(abs(func_p) > 10^-10) == 0,1,'last'), find(cumsum(abs(func_m) > 10^-10) == 0,1,'last'));
            right_index = max(numel(func_p) - find(cumsum(flipud(abs(func_p(:))) > 10^-10) == 0,1,'last'), numel(func_p) - find(cumsum(flipud(abs(func_m(:))) > 10^-10) == 0,1,'last'));
            
            
           
            t = linspace(obj.Omega.Points(left_index)-obj.Shift_Baths(alpha)*obj.Mu(alpha), obj.Omega.Points(right_index)-obj.Shift_Baths(alpha)*obj.Mu(alpha), 4000);
            
            func_p_equi = interp1(obj.Omega.Points-obj.Shift_Baths(alpha)*obj.Mu(alpha), func_p, t);
            func_m_equi = interp1(obj.Omega.Points-obj.Shift_Baths(alpha)*obj.Mu(alpha), func_m, t);
            
            pv_int_p = cauchy(t, func_p_equi, 10^-4);
            pv_int_m = cauchy(t, func_m_equi,10^-4);
            
            
            
            pv_int_p(isnan(pv_int_p)|isinf(pv_int_p)) = 0;
            pv_int_m(isnan(pv_int_m)|isinf(pv_int_m)) = 0;
            C = zeros(size(energies));
            C(dag_mat == 1) = spline(t, pv_int_p, energies(dag_mat==1));
            C(dag_mat == -1) = spline(t, pv_int_m, -energies(dag_mat==-1));
            
            
            
        end
        
        
        
        function [C, Bat, precalc_warning] = f_sigma(obj, omega, mu_alpha, beta_alpha, alpha, c_spin, dag)
            
            %TODO: introduce a general exception to trigger mes_flag
            mes_flag = false;
            
            
            % gamma = i/2 sum    GR-GA + sGK (dag\omega)   |_ \pi(dag)(lk)    ...   exchange l and k if dag is negative
            % in equilibrium K is given by (GR-GA) (1-2fermi)
            
            % sum = \sum_\alpha, kl   V_\alpha,k\mu^s V_\alpha,l\kappa^\overline(s)
            
            [n1, n2, npoints] = size(obj.GF{alpha,1});
            C = zeros(n1,n2,numel(omega));
            
            %shift baths works the following way:
            % 1) either: [omega, GF(omega) * fermi(omega, mu)] -> [omega, GF(omega+shift) * fermi(omega,mu)]
            % 2) or:     [omega, GF(omega) * fermi(omega, mu)] -> [omega - shift, GF(omega) * fermi(omega - shift,mu)]
            
            
           points = obj.Omega.Points - obj.Shift_Baths(alpha) * mu_alpha;
           %NOTE: due to the Fermi function just one half side of the integral gives nonzero values
           %(f^\overline{s}(\omega) > 10^-16!  ==>  log(10^16-1)/beta + dag * (points - mu) > 0
           L = log(10^16-1)/beta_alpha + dag*(points-mu_alpha)> 0;
           L = true(size(obj.Omega.Points));
           points = points(L);
           npoints = numel(points);
           
            % in bath equilibrium!
            GFAK =  dag/2*(obj.GF{alpha, 1}(:,:,L) - conj(permute(obj.GF{alpha,1}(:,:,L), [2,1,3]))) .* ...
                repmat(permute(1+dag-2*obj.Fermi(points, mu_alpha, beta_alpha, 1), [2,3,1]), n1, n2);
            %C = permute(interp1(obj.Omega.Points , permute(GFAK, [3,1,2]), dag*omega), [2,3,1]);
            for k = 1:n1
                for l = 1:n2
                    for en_in = 1: numel(omega)
                        
                        going_left = 7;
                        
                        left_index = find( points < dag * omega(en_in),1, 'last');
                        if isempty(left_index), left_index = 0; end
                        
                        [L_lower, L_upper, L_middle] = deal(true(1,npoints));
                        if left_index <= going_left +2
                           %  pole out of grid, to the left
                           L_lower = ~L_lower;
                           L_middle = ~L_middle;
                           treat_pole = false;
                        elseif left_index >= npoints-going_left-3
                            % pole outside on the right side
                           
                           L_middle = ~L_middle;
                           L_upper = ~L_upper;
                           treat_pole = false;
                        else
                            treat_pole = true;
                            L_lower(left_index-going_left:end) = false;
                            L_upper(1:left_index+going_left+1) = false;
                            L_middle(left_index-going_left-1:left_index+going_left+2) = false; L_middle = ~L_middle;
                        end
                        
                        
                        if treat_pole
                            polynom = spline(points(L_middle), squeeze(GFAK(k,l,L_middle)));
                            
                            
                            delta_x = dag*omega(en_in) - points(left_index - going_left-1); % for symmetric interval
                            pole_int = integral(@(x) (ppval(polynom,dag*(omega(en_in) - x)) - ppval(polynom, dag*(omega(en_in) + x)))./x, 0, delta_x) + ...
                                integral(@(x) ppval(polynom,x)./(omega(en_in) - dag*x), dag*omega(en_in) + delta_x, points(left_index+going_left+2));
                        else
                            pole_int = 0;
                        end
                        
                        C(k,l,en_in) = -1/pi * (trapz(points(L_lower)- 0*obj.Shift_Baths(alpha)*mu_alpha, squeeze(GFAK(k,l,L_lower))./(omega(en_in) - dag*points(L_lower))) + ...
                            trapz(points(L_upper), squeeze(GFAK(k,l,L_upper))./(omega(en_in) - dag*points(L_upper))) + pole_int);
                        % NOTE: Test if limits match
                        % [min(obj.Omega.Points(L_lower)), max(obj.Omega.Points(L_lower)), dag*omega(en_in)-delta_x, dag*omega(en_in) + delta_x, dag*omega(en_in)+delta_x, obj.Omega.Points(left_index+going_left+2), min(obj.Omega.Points(L_upper)), max(obj.Omega.Points(L_upper))]
                        
                        
                        % TEST
                        %poly = spline(points-obj.Shift_Baths(alpha)*mu_alpha, squeeze(GFAK(k,l,:))./(omega(en_in) - dag*points));
                        
                        %Test = -1/pi * integral(@(x) ppval(poly, x), -30,30);
                        
                        
                        % add integration of limits to -inf +inf by 
                        % Strategy: take last grid points and fit a 1/x function which can be integrated
                        % analytically.
                        L_left_boundary = 1:10;
                        L_right_boundary = numel(points) + (-9:0);
                        left_integrand =  squeeze(GFAK(k,l,L_left_boundary))./(omega(en_in) - dag*points(L_left_boundary));
                        right_integrand = squeeze(GFAK(k,l,L_right_boundary))./(omega(en_in) - dag*points(L_right_boundary));
                        
                      
                        left = true;
                        [left_integ,  left_parameters] = f_boundary_integral(points(L_left_boundary), left_integrand, left, mes_flag);
                        [right_integ, right_parameters] = f_boundary_integral(points(L_right_boundary), right_integrand, ~left, mes_flag);
                        
                        %plot to check convergence
%                         subplot(2,1,1)   
%                         plot([points(1)-1; points(L_left_boundary)], imag([0;left_integrand]), 'x')
%                         hold on
%                         plot([points(1)-1; points(L_left_boundary)], imag(1./([points(1)-1; points(L_left_boundary)] - left_parameters(2)).^2*left_parameters(1)))
%                         subplot(2,1,2)
%                         plot([points(L_left_boundary); points(end)+1], imag([left_integrand; 0]), '-x')
%                         hold on
%                         plot([points(L_left_boundary); points(end)+1], imag(1./([points(L_left_boundary); points(end)+1] - right_parameters(2)).^2*right_parameters(1)))


                        
                       
                        %NOTE: if 
                        %comparison plot
%                         t = linspace(-100,points(100), 100);
%                         
%                         plot(t, left_parameters(1)./t + left_parameters(2) * t, 'r-')
%                         hold on
%                         plot(points(1:100), squeeze(GFAK(k,l,1:100))./(omega(en_in) - dag*points(1:100)), 'xb')
%                         
%                         %right side
%                         t = linspace(points(end-100), 100,100);
%                         plot(t,right_parameters(1)./t + left_parameters(2) * t, 'b')
%                         hold on
%                         plot(points(end-100:end), squeeze(GFAK(k,l,end-100:end))./(omega(en_in) - dag*points(end-100:end)),'xr')
%                         
%                         plot(points, imag(squeeze(GFAK(k,l,:))), 'k.')

                         C(k,l,en_in) = C(k,l,en_in)  - 1/pi * (left_integ + right_integ); 
                        
                        
                    end
                end
            end
            %transpose C if dag is negative
            if dag == -1
                C = permute(C, [2,1,3]);
            end
            
            %the question is, if we should use precalculation? - maybe not...
            Bat = obj;
            precalc_warning = false;
            
              
            
        end
        
        
        function [C, obj, warn_precalc] = f_lookup_sigma(obj, energy, mu_alpha, beta_alpha, alpha, c_spin, dag)
            
            % calculates the integral:
            % int i(G^R(om) - G^A(om)) * n_alpha^s(om) /(-energy + s*om + i0^+)
            % for all band orbitals
            
            %NOTE: spinless case:
            
            if ~isempty(c_spin)
                if obj.N_spin == 1
                    warning('use [] for c_spin in spinless systems! or rethink bath definition of N_spin')
                end
            else
                % NOTE: default value for spinless systems!
                c_spin = -1;
            end
            
            warn_precalc = 0;
            %alpha and c_spin just transfer to GF - whether it's different or not
            
            % how to store entries?
            % energy ... real number
            % mu_alpha ... real number
            % beta_alpha ... real number, may be more discrete, may be a fixed array
            % alpha ... discrete number, is known before
            % c_spin ... discrete number (2 variables)
            % dag ... discrete number (2 variables)
            
            col_index = obj.g_column_index(mu_alpha, beta_alpha, alpha, c_spin, dag);
            %NOTE: col_index is a vector only if bath has orbitals (size of GF)
            row_index = obj.g_row_index(energy);
            %NOTE: implement cutoff due to limited band
            
            %for sigma there is no bandwidth as it goes with 1/x
            % Fermi function helps to distinguish
            % \Theta( s(\mu - energy) + log(10^16-1)/beta)
            
            L_ener = true(size(energy)); % dag*(energy - mu_alpha) + log(10^16-1)/beta_alpha > 0;
            
            
            %L_ener = (energy > obj.Omega.Limits(1)) & (energy < obj.Omega.Limits(2));
            %TODO: Why are those entries zero??? Poles outside of mesh where GFs are defined...
            % If GF are shifted, this becomes kind of relevant. Grid size should be big enough to have bandwidth
            % inside!!! -> so poles outside lead to zeros.
            
            %C = zeros(numel(L_ener), numel(col_index));
            n_orb = obj.Numel_Orbitals(alpha) * (2-obj.Spin_symmetric(alpha)); %sqrt(obj.Numel_Orbitals(alpha,1));
            C = zeros(n_orb, n_orb, numel(row_index));
            
            if any(isnan(col_index)) || any(isnan(row_index))
                % TODO: function for interpolate and update ranges (Beta_Range, En_Range, Mu_Range)
                warning('input values not preinitialized')
                warn_precalc = 1;
                
                C(:,:,L_ener) = obj.f_sigma(energy(L_ener), mu_alpha, beta_alpha, alpha, c_spin, dag);
                
            else
                
                %check if bath correlation function has been calculated:
                %if all(all(obj.Bath_Calc(row_index(L_ener), col_index) == true))
                %TODO: split, those which are already calculated shall be taken, the other shall be calculated
                % Check energies which are close enough
                
                
                % NOTE: define those values which shall be interpolated!!!
                L_alr_calc = full(all(obj.Bath_Calc(row_index,col_index) == 1,2)) & L_ener;
                
                En_range_calculated = obj.En_Range(all(obj.Bath_Calc(:,col_index) == 1,2));
                
                L_approx = L_ener & ~L_alr_calc;
                
                if sum(L_approx) > 0 && numel(En_range_calculated) > 10 && min(energy(L_ener)) > min(En_range_calculated) && max(energy(L_ener)) < max(En_range_calculated)
                    En_handle = griddedInterpolant(En_range_calculated, En_range_calculated, 'previous');
                    En_prev = En_handle(energy(L_ener));
                    En_handle.Method = 'next';
                    
                    En_next = En_handle(energy(L_ener));
                    En_handle.Method = 'nearest';
                    %[En_prev, energy(L_ener), En_next]  % all three are equal if that
                    
                    %L_approx(L_ener) = (En_next - En_prev) < obj.Omega.Delta/500   & L_approx(L_ener);
                    L_approx(L_ener) = (En_next - En_prev) < 10^-7   & L_approx(L_ener);
                    if obj.Omega.Uniform == false && sum(L_approx) > 0
                        %warning('Grid not uniform, Delta not well defined!')
                        
                    end
                    
                    %TODO: Delta may not be well defined for not uniform grids!
                else
                    L_approx = L_ener*0;
                end
                %En_nearest = En_handle(energy(L_ener));
                % [En_prev, En_next, abs(energy(L_ener) - En_nearest), (En_next - En_prev) < obj.Omega.Delta, energy(L_ener), full(all(obj.Bath_Calc(row_index,col_index) == 1,2))]
                
                %TODO: Strategy: identify intervals, which are so close, that delta is smaller than Grid Delta
                %and use linear interpolation for those.
                
                %disable interpolation:
                L_approx = L_ener*0;
              
                
                % NOTE: define which values are calculated
                L_calc = L_ener & ~ L_alr_calc & ~L_approx;

                               
                
                % NOTE: already calculated
                C_temp = full(obj.Bath_Corr(row_index(L_alr_calc), col_index));
                
                C(:,:, L_alr_calc) = reshape(C_temp.', n_orb, n_orb, sum(L_alr_calc));
                
                % linear approximation for each orbital
                if any(L_approx)
                    disp([num2str(sum(L_approx)), ' values of bath correlation function approximated. Delta: ', num2str(obj.Omega.Delta)])
                    C_temp = zeros(sum(L_approx), numel(col_index));
                    for k = 1: numel(col_index)
                        
                        C_temp(:,k) = interp1(En_range_calculated, full(obj.Bath_Corr(obj.Bath_Calc(:,col_index(k)) == 1, col_index(k))), energy(L_approx), 'linear');
                        
                    end
                    
                    C(:,:, L_approx) = reshape(C_temp.', n_orb, n_orb, sum(L_approx));
                end
                
                %NOTE: old version
                %                     C_temp = full(obj.Bath_Corr(row_index(L_ener), col_index));
                %                     n_orb = sqrt(obj.Numel_Orbitals(alpha,obj.Spin_fun(c_spin)));
                %
                %                     C(:,:, L_ener) = reshape(C_temp.', n_orb, n_orb, sum(L_ener));
                
                
                %else
                %TODO: calculate the integral
                
              
                
                
                C(:,:,L_calc) = obj.f_sigma(energy(L_calc), mu_alpha, beta_alpha, alpha, c_spin, dag);
                       
                obj.Bath_Corr(row_index, col_index) = reshape(C,obj.Numel_Orbitals(alpha,1)^2, numel(energy)).';
                obj.Bath_Calc(row_index, col_index) = true;
                
                %end
            end
            
            % use Tab
            % C = obj.Categ.f_access('Corr_Bath', 'EMBASD', energy, mu_alpha, beta_alpha, alpha, c_spin, dag);
            % if isempty(C)
            % alternative would be preallocation - squared list...
            % other alternative would be sparse matrix, preindexed by list of energies etc.
            % corr_bath = obj.f_calculate(energy, mu_alpha, beta_alpha, alpha, c_spin, dag);
            % obj.Categ.f_add('EMBASDC', energy, mu_alpha, beta_alpha, alpha, c_spin, dag, corr_bath);
            % C = corr_bath;
            % end
        end
        
        function obj = f_refine_GF(obj, new_grid)
            if ~isempty(obj.GF)
                
                % for all baths all spins, all Bath orbitals_row all Bath orbitals_column
                
                if ~isa(new_grid, 'Grid')
                    new_grid = Grid(new_grid);
                end
                
                for k = 1:obj.Numel_Baths
                    
                    
                    obj.GF{k,1} = spline(obj.Omega.Points, obj.GF{k,1}, new_grid.Points);
                    
                    
                    
                end
                obj.Omega = new_grid;
                
            else
                error('no Green function defined!')
            end
            
        end
        
        function [C, obj, warn_precalc] = f_lookup(obj, energy, mu_alpha, beta_alpha, alpha, c_spin, dag)
            
            % calculates the integral:
            % int i(G^R(om) - G^A(om)) * n_alpha^s(om) /(-energy + s*om + i0^+)
            % for all band orbitals
            
            %NOTE: spinless case:
            
            if ~isempty(c_spin)
                if obj.Numel_Spins > 1
                    warning('use [] for c_spin in spinless systems! or rethink bath definition of N_spin')
                end
            else
                % NOTE: default value for spinless systems!
                c_spin = -1;
            end
            
            warn_precalc = 0;
            %alpha and c_spin just transfer to GF - whether it's different or not
            
            % how to store entries?
            % energy ... real number
            % mu_alpha ... real number
            % beta_alpha ... real number, may be more discrete, may be a fixed array
            % alpha ... discrete number, is known before
            % c_spin ... discrete number (2 variables)
            % dag ... discrete number (2 variables)
            
            col_index = obj.g_column_index(mu_alpha, beta_alpha, alpha, c_spin, dag);
            %NOTE: col_index is a vector only if bath has orbitals (size of GF)
            row_index = obj.g_row_index(energy);
            %NOTE: implement cutoff due to limited band
            
            L_ener = (energy > obj.Omega.Limits(1)) & (energy < obj.Omega.Limits(2));
            %TODO: Why are those entries zero??? Poles outside of mesh where GFs are defined...
            % If GF are shifted, this becomes kind of relevant. Grid size should be big enough to have bandwidth
            % inside!!! -> so poles outside lead to zeros.
            
            %C = zeros(numel(L_ener), numel(col_index));
            n_orb = sqrt(obj.Numel_Orbitals(alpha,obj.Spin_fun(c_spin)));
            C = zeros(n_orb, n_orb, numel(row_index));
            
            if any(isnan(col_index)) || any(isnan(row_index))
                % TODO: function for interpolate and update ranges (Beta_Range, En_Range, Mu_Range)
                warning('input values not preinitialized')
                warn_precalc = 1;
                if size(obj.GF{alpha,obj.Spin_fun(c_spin)},3) == 1
                    if ~isdeployed()
                        % disp('Flat band approximation')
                    end
                    C(:,:,L_ener) = obj.f_flat_band_calculate(energy(L_ener), mu_alpha, beta_alpha, alpha, c_spin, dag);
                else
                    switch obj.Integrator
                        case 'Fast'
                            C(:,:,L_ener) = obj.f_calculate_fast(energy(L_ener), mu_alpha, beta_alpha, alpha, c_spin, dag);
                        case 'Matlab'
                            C(:,:,L_ener) = obj.f_calculate(energy(L_ener), mu_alpha, beta_alpha, alpha, c_spin, dag);
                        case 'Detail'
                            C(:,:,L_ener) = obj.f_calculate_detail(energy(L_ener), mu_alpha, beta_alpha, alpha, c_spin, dag);
                        case 'Fast_Detail'
                            C(:,:,L_ener) = obj.f_calculate_detail_fast(energy(L_ener), mu_alpha, beta_alpha, alpha, c_spin, dag);
                        case 'Constant'
                            C(:,:,L_ener) = obj.f_flat_band_calculate(energy(L_ener), mu_alpha, beta_alpha, alpha, c_spin, dag);
                            
                    end
                end
            else
                
                %check if bath correlation function has been calculated:
                %if all(all(obj.Bath_Calc(row_index(L_ener), col_index) == true))
                %TODO: split, those which are already calculated shall be taken, the other shall be calculated
                % Check energies which are close enough
                
                
                % NOTE: define those values which shall be interpolated!!!
                L_alr_calc = full(all(obj.Bath_Calc(row_index,col_index) == 1,2)) & L_ener;
                
                En_range_calculated = obj.En_Range(all(obj.Bath_Calc(:,col_index) == 1,2));
                
                L_approx = L_ener & ~L_alr_calc;
                
                if sum(L_approx) > 0 && numel(En_range_calculated) > 10 && min(energy(L_ener)) > min(En_range_calculated) && max(energy(L_ener)) < max(En_range_calculated)
                    En_handle = griddedInterpolant(En_range_calculated, En_range_calculated, 'previous');
                    En_prev = En_handle(energy(L_ener));
                    En_handle.Method = 'next';
                    
                    En_next = En_handle(energy(L_ener));
                    En_handle.Method = 'nearest';
                    %[En_prev, energy(L_ener), En_next]  % all three are equal if that
                    
                    L_approx(L_ener) = (En_next - En_prev) < obj.Omega.Delta   & L_approx(L_ener);
                    if obj.Omega.Uniform == false && sum(L_approx) > 0
                        %warning('Grid not uniform, Delta not well defined!')
                        
                    end
                    
                    %TODO: Delta may not be well defined for not uniform grids!
                else
                    L_approx = L_ener*0;
                end
                %En_nearest = En_handle(energy(L_ener));
                % [En_prev, En_next, abs(energy(L_ener) - En_nearest), (En_next - En_prev) < obj.Omega.Delta, energy(L_ener), full(all(obj.Bath_Calc(row_index,col_index) == 1,2))]
                
                %TODO: Strategy: identify intervals, which are so close, that delta is smaller than Grid Delta
                %and use linear interpolation for those.
                
                %disable interpolation:
                L_approx = L_ener*0;
              
                
                % NOTE: define which values are calculated
                L_calc = L_ener & ~ L_alr_calc & ~L_approx;

                               
                
                % NOTE: already calculated
                C_temp = full(obj.Bath_Corr(row_index(L_alr_calc), col_index));
                n_orb = sqrt(obj.Numel_Orbitals(alpha,obj.Spin_fun(c_spin)));
                
                C(:,:, L_alr_calc) = reshape(C_temp.', n_orb, n_orb, sum(L_alr_calc));
                
                % linear approximation for each orbital
                if any(L_approx)
                    disp([num2str(sum(L_approx)), ' values of bath correlation function approximated. Delta: ', num2str(obj.Omega.Delta)])
                    C_temp = zeros(sum(L_approx), numel(col_index));
                    for k = 1: numel(col_index)
                        
                        C_temp(:,k) = interp1(En_range_calculated, full(obj.Bath_Corr(obj.Bath_Calc(:,col_index(k)) == 1, col_index(k))), energy(L_approx), 'linear');
                        
                    end
                    
                    C(:,:, L_approx) = reshape(C_temp.', n_orb, n_orb, sum(L_approx));
                end
                
                %NOTE: old version
                %                     C_temp = full(obj.Bath_Corr(row_index(L_ener), col_index));
                %                     n_orb = sqrt(obj.Numel_Orbitals(alpha,obj.Spin_fun(c_spin)));
                %
                %                     C(:,:, L_ener) = reshape(C_temp.', n_orb, n_orb, sum(L_ener));
                
                
                %else
                %TODO: calculate the integral
                
                %distinguish between integral and flat band approximation
                
                if size(obj.GF{alpha,obj.Spin_fun(c_spin)},3) == 1
                    if ~isdeployed()
                        %disp('Flat band approximation')
                    end
                    C(:,:,L_calc) = obj.f_flat_band_calculate(energy(L_calc), mu_alpha, beta_alpha, alpha, c_spin, dag);
                else
                    switch obj.Integrator
                        case 'Fast'
                            C(:,:,L_calc) = obj.f_calculate_fast(energy(L_calc), mu_alpha, beta_alpha, alpha, c_spin, dag);
                        case 'Matlab'
                            C(:,:,L_calc) = obj.f_calculate(energy(L_calc), mu_alpha, beta_alpha, alpha, c_spin, dag);
                        case 'Detail'
                            C(:,:,L_calc) = obj.f_calculate_detail(energy(L_calc), mu_alpha, beta_alpha, alpha, c_spin, dag);
                        case 'Fast_Detail'
                            C(:,:,L_calc) = obj.f_calculate_detail_fast(energy(L_calc), mu_alpha, beta_alpha, alpha, c_spin, dag);
                        case 'Constant'
                            C(:,:,L_calc) = obj.f_flat_band_calculate(energy(L_calc), mu_alpha, beta_alpha, alpha, c_spin, dag);
                            
                            
                    end
                    
                    
                end
                obj.Bath_Corr(row_index, col_index) = reshape(C,obj.Numel_Orbitals(alpha,obj.Spin_fun(c_spin)), numel(energy)).';
                obj.Bath_Calc(row_index, col_index) = true;
                
                %end
            end
            
            % use Tab
            % C = obj.Categ.f_access('Corr_Bath', 'EMBASD', energy, mu_alpha, beta_alpha, alpha, c_spin, dag);
            % if isempty(C)
            % alternative would be preallocation - squared list...
            % other alternative would be sparse matrix, preindexed by list of energies etc.
            % corr_bath = obj.f_calculate(energy, mu_alpha, beta_alpha, alpha, c_spin, dag);
            % obj.Categ.f_add('EMBASDC', energy, mu_alpha, beta_alpha, alpha, c_spin, dag, corr_bath);
            % C = corr_bath;
            % end
        end
        
        function [C, func] = f_calculate(obj, energ_inp, mu_alpha, beta_alpha, alpha, c_spin, dag)
            C = zeros(size(obj.GF{alpha, 1},1), size(obj.GF{alpha, 1},2), numel(energ_inp));
            
            
            
            for row = 1:size(obj.GF{alpha, 1}, 1)
                for column = 1:size(obj.GF{alpha, 1}, 2)
                    if sum(abs(obj.GF{alpha,1}(row, column,:))) > 10^-10
                        
                        pp = spline(obj.Omega.Points, squeeze(obj.GF{alpha, 1}(row, column, :)) - conj(squeeze(obj.GF{alpha, 1}(column, row, :))));
                        
                        for en_ind = 1:numel(energ_inp)
                            energy = energ_inp(en_ind);
                            
                            func = @(x) - ppval(pp,x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(x,mu_alpha, beta_alpha, dag)./(2*pi *(dag * x + energy + 1i*obj.Omega.Imag_Infinitesimal));
                            
                            %C(-E)
                            % NOTE: Bath shift by mu_alpha optionally included
                            C(row, column, en_ind) = integral(@(x) - ppval(pp,x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(x,mu_alpha, beta_alpha, dag)./(2*pi *(dag * x + energy + 1i*obj.Omega.Imag_Infinitesimal)),...
                                obj.Omega.Limits(1), obj.Omega.Limits(2));
                            
                        end
                    else
                        C(row, column, :) = zeros(numel(energ_inp),1);
                    end
                end
            end
            
        end
        
        function C = f_calculate_detail(obj, energ_inp, mu_alpha, beta_alpha, alpha, c_spin, dag)
            C = zeros(size(obj.GF{alpha, 1},1), size(obj.GF{alpha, 1},2), numel(energ_inp));
            
            epsi = 10^-6;%obj.Omega.Imag_Infinitesimal;
            
            for row = 1:size(obj.GF{alpha, 1}, 1)
                for column = 1:size(obj.GF{alpha, 1}, 2)
                    if sum(abs(obj.GF{alpha,1}(row, column,:))) > 10^-10
                        
                        ZZ = zeros(size(energ_inp));
                        pp = spline(obj.Omega.Points, squeeze(obj.GF{alpha,1}(row, column, :)) - conj(squeeze(obj.GF{alpha, 1}(column, row, :))));
                        
                        for en_ind = 1:numel(energ_inp)
                            energy = energ_inp(en_ind);
                            
                            %C(-E)
                            
                            % NOTE: Bath shift by mu_alpha optionally included
                            ZZ(en_ind) = integral(@(x) - ppval(pp,x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(x,mu_alpha, beta_alpha, dag)./(2*pi *(dag * x + energy)),...
                                obj.Omega.Limits(1), -dag * energy - epsi) + ...
                                integral(@(x) - ppval(pp,x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(x,mu_alpha, beta_alpha, dag)./(2*pi *(dag * x + energy)),...
                                -dag * energy + epsi, obj.Omega.Limits(2)) + ...
                                integral(@(x) ((- ppval(pp,-dag*energy + x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(-dag*energy + x,mu_alpha, beta_alpha, dag))   - ...
                                (- ppval(pp,-dag*energy - x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(-dag*energy - x,mu_alpha, beta_alpha, dag))) ./(x*dag), 0, epsi);
                            
                            
                        end
                        
                        C(row,column,:) =  (ZZ -  1i/2 * -ppval(pp, -dag*energ_inp - obj.Shift_Baths(alpha)*mu_alpha).*obj.Fermi(-dag*energ_inp, mu_alpha, beta_alpha, dag)); %the i comes from performing the laplace integral at om = 0
                    else
                        C(row,column,:) = zeros(numel(energ_inp),1);
                    end
                end
            end
            
        end
        
        function C = f_calculate_detail_fast(obj, energ_inp, mu_alpha, beta_alpha, alpha, c_spin, dag)
            C = zeros(size(obj.GF{alpha, 1},1), size(obj.GF{alpha,1},2), numel(energ_inp));
            %aim is to do a rough approximation outside the function (epsilon environment 
            % the question is how big the relevant area is.
            %epsi = 10^-4;
            %error('does not work!!!')
            %TODO: comparison with detail integration shows severe problems!!!
            
            
            for row = 1:size(obj.GF{alpha, 1}, 1)
                for column = 1:size(obj.GF{alpha, 1}, 2)
                    if sum(abs(obj.GF{alpha,1}(row, column,:))) > 10^-10
                        ZZ = zeros(size(energ_inp));
                        pp = spline(obj.Omega.Points,  squeeze(obj.GF{alpha, 1}(row, column, :)) - conj(squeeze(obj.GF{alpha, 1}(column, row, :))));
                        
                        for en_ind = 1:numel(energ_inp)
                            energy = energ_inp(en_ind);
                            
                            pole = -dag*energy;
                            %TODO: Problem with eps-environment (spline area around!!!), eps should be around 10^-3
                            int_node_length = 3;
                            index = obj.Omega.f_find_index(pole) + [-int_node_length,int_node_length];
                            if any(index < 3) || any(index > (obj.Omega.N-2))
                                warning('pole out of grid - return zero')
                                if ~isdeployed()
                                    %disp('pole out of grid')
                                end
                            else
                                epsi = pole - obj.Omega.Points(index(1)); %symmetric distance
                                %case where pole (o) is closer to the first grid index(1), necessary to integrate from | to second X
                                %--------x--------X_____o__i__|=====X--------x--------x
                                %case where polo (o) is closer to the second grid index(2), do index(2) + 1 and
                                %integrate from | to new X
                                %--------x--------X________i_o______x___|====X
                                % ___ symmetric interval
                                if pole > obj.Omega.Points(index(1)+int_node_length)
                                    index(2) = index(2)+1;
                                end
                                
                                
                                first_range = 1:index(1);
                                x = obj.Omega.Points(first_range);
                                x_shift = obj.Omega.Points(first_range) - obj.Shift_Baths(alpha)*mu_alpha;
                                
                                ZZ(en_ind) = trapz(x_shift, -(squeeze(obj.GF{alpha, 1}(row, column, first_range)) - conj(squeeze(obj.GF{alpha, 1}(column, row, first_range))))...
                                    .* obj.Fermi(x,mu_alpha, beta_alpha, dag)./(2*pi *(dag * x + energy )) );
                                
                                second_range = index(2):obj.Omega.N;
                                x2 = obj.Omega.Points(second_range);
                                x2_shift = obj.Omega.Points(second_range) - obj.Shift_Baths(alpha)*mu_alpha;
                                
                                ZZ(en_ind) = ZZ(en_ind) + trapz(x2_shift, -(squeeze(obj.GF{alpha, 1}(row, column, second_range)) - conj(squeeze(obj.GF{alpha, 1}(column, row, second_range))))...
                                    .* obj.Fermi(x2,mu_alpha, beta_alpha, dag)./(2*pi *(dag * x2 + energy )) );
                                
                                middle_range = (index(1)) : (index(2));
                                x_middle = obj.Omega.Points(middle_range);
                                %pp = spline(x_middle,  squeeze(obj.GF{alpha, obj.Spin_fun(c_spin)}(row, column, middle_range)) - conj(squeeze(obj.GF{alpha, obj.Spin_fun(c_spin)}(column, row, middle_range))));
                                
                                %C(-E)
                                
                                % NOTE: Bath shift by mu_alpha optionally included
                                %TODO: first two integrals use trapz!!!
                                %ZZ(en_ind) = ZZ(en_ind) + ...
                                %   integral(@(x) ((- ppval(pp,-dag*energy + x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(-dag*energy + x,mu_alpha, beta_alpha, dag))   - ...
                                %  (- ppval(pp,-dag*energy - x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(-dag*energy - x,mu_alpha, beta_alpha, dag))) ./(x*dag), 0, epsi);
                                
                                ZZ(en_ind) = ZZ(en_ind) + ...
                                    integral(@(x) ((- ppval(pp,-dag*energy + x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(-dag*energy + x,mu_alpha, beta_alpha, dag))   - ...
                                    (- ppval(pp,-dag*energy - x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(-dag*energy - x,mu_alpha, beta_alpha, dag))) ./(x*dag), 0, epsi);
                                
                                %integrate rest between symmetric interval and right trapz region:
                                ZZ(en_ind) = ZZ(en_ind) + ...
                                    integral(@(x)- ppval(pp,x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(x,mu_alpha, beta_alpha, dag)./(2*pi *(dag * x + energy)),...
                                    pole + epsi, obj.Omega.Points(index(2)));
                                
                                
                                
                                %diff([x_shift(1),x_shift(end), pole-epsi-obj.Shift_Baths(alpha)*mu_alpha, pole+epsi-obj.Shift_Baths(alpha)*mu_alpha, pole+epsi-obj.Shift_Baths(alpha)*mu_alpha, obj.Omega.Points(index(2))-obj.Shift_Baths(alpha)*mu_alpha, min(x2_shift), max(x2_shift)])
                                
                            end
                            
                        end
                        %NOTE: estimate when functions are out of GF range to contribute.
                        %NOTE: find a supplement of spline, since it costs a lot
                        %NOTE:
                        C(row,column,:) =  (ZZ -  1i/2 * -ppval(pp, -dag*energ_inp - obj.Shift_Baths(alpha)*mu_alpha).*obj.Fermi(-dag*energ_inp, mu_alpha, beta_alpha, dag)); %the i comes from performing the laplace integral at om = 0
                    else
                        C(row,column,:) = zeros(numel(energ_inp),1);
                    end
                end
            end
            
        end
        
        function C = f_calculate_detail_fast_2(obj, energ_inp, mu_alpha, beta_alpha, alpha, c_spin, dag)
            warning('not so exact')
            C = zeros(size(obj.GF{alpha, 1},1), size(obj.GF{alpha, 1},2), numel(energ_inp));
            %aim is to do a rough approximation outside the function (epsilon environment 
            % the question is how big the relevant area is.
            %epsi = 10^-4;
            
            for row = 1:size(obj.GF{alpha, 1}, 1)
                for column = 1:size(obj.GF{alpha, 1}, 2)
                    ZZ = zeros(size(energ_inp));
                    pp = spline(obj.Omega.Points,  squeeze(obj.GF{alpha, 1}(row, column, :)) - conj(squeeze(obj.GF{alpha, 1}(column, row, :))));
                        
                    for en_ind = 1:numel(energ_inp)
                        energy = energ_inp(en_ind);
                        
                        pole = -dag*energy;
                        %TODO: Problem with eps-environment (spline area around!!!), eps should be around 10^-3
                        index = obj.Omega.f_find_index(pole) + [-1,1];
                        if any(index < 3) || any(index > (obj.Omega.N-2))
                            if ~isdeployed()
                                disp('pole out of grid')
                                %NOTE: is it is out of grid, one should exclude it before the for loop?
                            end
                        else
                            epsi = diff(obj.Omega.Points(index));
                            
                            first_range = 1:index(1);
                            x = obj.Omega.Points(first_range);
                            x_shift = obj.Omega.Points(first_range) - obj.Shift_Baths(alpha)*mu_alpha;
                            
                            ZZ(en_ind) = trapz(x_shift, -(squeeze(obj.GF{alpha, 1}(row, column, first_range)) - conj(squeeze(obj.GF{alpha, 1}(column, row, first_range))))...
                                .* obj.Fermi(x,mu_alpha, beta_alpha, dag)./(2*pi *(dag * x + energy )) );
                            
                            second_range = index(2):obj.Omega.N;
                            x2 = obj.Omega.Points(second_range);
                            x2_shift = obj.Omega.Points(second_range) - obj.Shift_Baths(alpha)*mu_alpha;
                            
                            ZZ(en_ind) = ZZ(en_ind) + trapz(x2_shift, -(squeeze(obj.GF{alpha, 1}(row, column, second_range)) - conj(squeeze(obj.GF{alpha, 1}(column, row, second_range))))...
                                .* obj.Fermi(x2,mu_alpha, beta_alpha, dag)./(2*pi *(dag * x2 + energy )) );
                            
                            middle_range = (index(1)) : (index(2));
                            x_middle = obj.Omega.Points(middle_range);
                            %pp = spline(x_middle,  squeeze(obj.GF{alpha, obj.Spin_fun(c_spin)}(row, column, middle_range)) - conj(squeeze(obj.GF{alpha, obj.Spin_fun(c_spin)}(column, row, middle_range))));
                            
                            %C(-E)
                            
                            % NOTE: Bath shift by mu_alpha optionally included
                            %TODO: first two integrals use trapz!!! this one consists of three points - use
                            %simpson!!!
                            
                            %tt = linspace(min(x_middle), max(x_middle), 100);
                            tt = linspace(0, epsi, 100);
                            func =  @(x) ((- ppval(pp,-dag*energy + x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(-dag*energy + x,mu_alpha, beta_alpha, dag))   - ...
                                (- ppval(pp,-dag*energy - x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(-dag*energy - x,mu_alpha, beta_alpha, dag))) ./(x*dag);
                            plot(tt, real(func(tt)), tt, imag(func(tt)))
                            
                            
                            ZZ(en_ind) = ZZ(en_ind) + ...
                                integral(@(x) ((- ppval(pp,-dag*energy + x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(-dag*energy + x,mu_alpha, beta_alpha, dag))   - ...
                                (- ppval(pp,-dag*energy - x-obj.Shift_Baths(alpha)*mu_alpha) .* obj.Fermi(-dag*energy - x,mu_alpha, beta_alpha, dag))) ./(x*dag), 0, epsi);
                        end
                        
                    end
                    %NOTE: estimate when functions are out of GF range to contribute.
                    %NOTE: find a supplement of spline, since it costs a lot
                    %NOTE: 
                    C(row,column,:) =  (ZZ -  1i/2 * -ppval(pp, -dag*energ_inp - obj.Shift_Baths(alpha)*mu_alpha).*obj.Fermi(-dag*energ_inp, mu_alpha, beta_alpha, dag)); %the i comes from performing the laplace integral at om = 0
                end
            end
            
        end
        
        function C = f_calculate_fast(obj, energ_inp, mu_alpha, beta_alpha, alpha, c_spin, dag)
            C = zeros(size(obj.GF{alpha, 1},1), size(obj.GF{alpha, 1},2));
            warning('not very accurate')
            
            
            for row = 1:size(obj.GF{alpha, 1}, 1)
                for column = 1:size(obj.GF{alpha, 1}, 2)
                    for en_ind = 1:numel(energ_inp)
                        energy = energ_inp(en_ind);
                        %shift of green's function equal to -shift of Omega.Points
                        x = (obj.Omega.Points + obj.Shift_Baths(alpha)*mu_alpha);
                        C(row, column, en_ind) = trapz(x, -(squeeze(obj.GF{alpha, 1}(row, column,:)) - conj(squeeze(obj.GF{alpha, 1}(column, row,:)) )) ...
                            .* obj.Fermi(x,mu_alpha, beta_alpha, dag)./(2*pi *(dag * x + energy + 1i*obj.Omega.Imag_Infinitesimal)) );
                        
                        
                    end
                end
            end
            
        end
        
        function C = f_flat_band_calculate(obj, energ_inp, mu_alpha, beta_alpha, alpha, c_spin, dag, eps_0)
            if nargin < 8 || isempty(eps_0), eps_0 = 0; end
            C = zeros(size(obj.GF{alpha,1},1), size(obj.GF{alpha, 1},2));
            
             for en_ind = 1:numel(energ_inp)
                energy = energ_inp(en_ind);
                
                for row = 1:size(obj.GF{alpha, 1}, 1)
                    for column = 1:size(obj.GF{alpha, 1}, 2)
                        
                        

                        %TODO: add offset (shift_bath) with Heaviside functions!
                        %NOTE: implemented by full G^ret function, otherwise just simple
                        
                        % defined as C(-energy)!!!
                        % NOTE: Bath shift by mu_alpha optionally included
                        constant = 1i * (obj.f_G_eval(eps_0-obj.Shift_Baths(alpha)*mu_alpha, alpha, c_spin, row, column) - conj(obj.f_G_eval(eps_0-obj.Shift_Baths(alpha)*mu_alpha, alpha, c_spin, column, row))) / (2*pi);
                        

                        
                        if abs(imag(constant) ) > 10^-10
                             warning('mixing of orbitals - check it') % real parts of Green functions should cancel
                        end
                        
                        %NOTE: if you use Fermi(-dag) then you have to give the imaginary part a minus sign too!
                        C(row, column, en_ind) = constant * ( pi * obj.Fermi(-dag*energy, mu_alpha, beta_alpha, dag) +...
                                                              1i * real(psiz(1/2+1i*beta_alpha*(energy + dag*mu_alpha)/(2*pi))));
                        % compare:
                        %ZZ = 2*pi*1./(pi*abs(obj.TB_hopping))* (1/2* 1./(exp(tau*BB_inp.*(-tau*EE_inp-MM_inp))+1) + ...
                        %                                    1i/(2*pi) * real(psiz(1/2 + 1i/(2*pi) * BB_inp.*(EE_inp + tau*MM_inp))));
                        
                        %NOTE: corresponds to F(-energy*dag, mu_alpha, beta_alpha)/2
                
      
                        
                    end
                end
            end
        
        
        end
        %just if Table is used...
        %function C = Corr_Bath(obj,varargin)
        %    C = obj.Categ.f_access('Corr_Bath', varargin{:});
        %end
        
        
        function obj =  f_preinitialize(obj, En, mu_vec, beta_vec, energy_cut)
            %NOTE: preinitialized for one setup
            beta_vec = unique(beta_vec(:));
            mu_vec = unique(mu_vec(:));
            
            Delta_En = obj.f_find_energy_delta(En, energy_cut);
            
            numel_orbitals = max(cellfun(@(x) size(x,1), obj.GF).^2); %takes the maximum of all orbitals to have enough space
            
            relevant_energy_states = numel(unique([-Delta_En, Delta_En]));
            relevant_Vb_dag_spin_states = numel(unique([-mu_vec,mu_vec])) * obj.Numel_Baths * numel_orbitals * 4 * numel(beta_vec);
            
            
            obj.Bath_Corr = sparse(relevant_energy_states, relevant_Vb_dag_spin_states);
            obj.Bath_Calc = sparse(relevant_energy_states, relevant_Vb_dag_spin_states); %stores a logical value if Calculation has been performed
            
            obj.En_Range = unique([-Delta_En, Delta_En]);
            obj.Mu_Range = mu_vec; %unique([-mu_vec,mu_vec]);
            obj.Beta_Range = beta_vec;
            
            
            
            
        end
        
        function row_index = g_row_index(obj,energy)
            [found, row_index] = (ismember(energy,obj.En_Range));
            if ~all(found)
                row_index = nan(size(found));
            end
        end
        
        function col_index = g_column_index(obj, mu_alpha, beta_alpha, alpha, c_spin, dag)
            %returns nan if index is not found
            
           
            
            
            mu_ind = find(mu_alpha == obj.Mu_Range);
            beta_ind = find(beta_alpha == obj.Beta_Range);
            
            spin_dag_ind = c_spin/2+1.5 + (dag/2+0.5) * 2; 
            % dag  spin value
            % -1    -1    1
            % -1    +1    2
            %  1    -1    3
            %  1    +1    4
            
            beta_max = numel(obj.Beta_Range);
            
            bath_max = obj.Numel_Baths;
           
            if isempty(mu_ind) || isempty(beta_ind)
                col_index = nan;
            else
            
            %         mu,    beta,    alpha,   spin_dag,   
                col_index =  ((((mu_ind - 1) * beta_max  + (beta_ind - 1)) * bath_max + (alpha - 1)) * 4   + (spin_dag_ind-1) )* max(cellfun(@(x) size(x,1), obj.GF).^2) + (1:(size(obj.GF{alpha,1},1) * size(obj.GF{alpha,1},2)));
            end
        end
        
        
        function f_find_index
            %return the parameters where everything is stored!
       
        end
        
        
        
        function [diff, ddd] = f_find_energy_delta(~, En, energy_cut)
            %TODO: check if all relevant energies are calculated!!!!
            % this function finds all relevant energy differences which are possible
            if nargin < 3 || isempty(energy_cut), energy_cut_usage = false; else, energy_cut_usage = true; end
            diff = [];
            if ismember('A', En.Structure)
                tag = 'NSA'; 
            else
                tag = 'NSE';
            end
            
            if En.Basis.N_spin == 1
                spinless = true;
            else
                spinless = false;
            end
            
            
            for N_el = 0:En.Basis.N_site*En.Basis.N_spin
                if spinless
                    Spin_range = 1;
                else
                    Spin_range = En.Basis.f_spin_range(N_el);
                end
                for Spin_ind = 1:numel(Spin_range)
                    if ~spinless
                        Spin = Spin_range(Spin_ind);
                    else
                        Spin = [];
                    end
 
                    
                    if energy_cut_usage
                        
                        A = unique(En.Energies(tag, N_el, Spin, {-inf, energy_cut}));
                    
                        up = unique(En.Energies('NS', N_el + 1, Spin + 1));
                        %up = [up,  unique(En.Energies('NSE', N_el - 1, Spin + 1, {-10000, energy_cut}))
                        down = unique(En.Energies('NS', N_el+1, Spin - 1));
                        up_ =  unique(En.Energies('NS', N_el - 1, Spin + 1));
                        down_ = unique(En.Energies('NS', N_el-1, Spin - 1));

                    else
                        A = unique(En.Energies('NS', N_el, Spin));
                    
                        up = unique(En.Energies('NS', N_el + 1, Spin + 1));
                        down = unique(En.Energies('NS', N_el+1, Spin - 1));
                        up_ = [];
                        down_ = [];
                    end
                    
                    [A1,A2] = meshgrid(A,[up; up_]);
                    [B1,B2] = meshgrid(A,[down; down_]);
                    
                    diff = [diff; abs(A1(:) - A2(:)); abs(B1(:) - B2(:))];
                    
                    
                    
                end
            end
            
            diff = unique(diff);
            
            
            ddd = 0;
            %[AA, BB] = meshgrid(unique(En.Energies));
            %ddd  = unique(abs(AA(:)-BB(:)));
            
            
            
            
        end
        
        
        function s_plot_energy_delta(obj, ddd)
            
            
            dif = diff(ddd);
            
            semilogy(ddd(1:end-1),dif)
        
        
        
            
        end
        
        function s_plot_bath_calc(obj)
           
            spy(obj.Bath_Calc)
        end
        
        function s_plot_GF(obj,bath, plot_spin, orbital_id)
            if nargin < 2 || isempty(bath), bath = 1; end
            if nargin < 3 || isempty(plot_spin), plot_spin = 1; end
            if nargin < 4 || isempty(orbital_id), orbital_id = 1; end
            %plots the Bath greens functions
            %{bath,spin}(orbital, orbital)
            obj.Numel_Orbitals, 
            obj.Numel_Baths
            obj.GF
            
          
            
            g_rows = numel(plot_spin) * 2; %times two for real and imaginary part
            g_columns = numel(orbital_id); 
            for k = 1:numel(plot_spin)
                 for l = 1:g_columns
                     
                     [ii,jj] = ind2sub([obj.Numel_Orbitals, obj.Numel_Orbitals], orbital_id(l));
                     %real
                     subplot(g_rows, g_columns,(l-1)*g_rows + 2*(k-1)  + 1)
                        hold on
                        plot(obj.Omega.Points, real(squeeze(obj.GF{bath, k}( ii, jj,:))))
                     %imag
                     subplot(g_rows, g_columns,(l-1)*g_rows + 2*(k-1)  + 2)
                        hold on
                        plot(obj.Omega.Points, imag(squeeze(obj.GF{bath, k}( ii, jj,:))))
                 end
            end
        end
        
        
        function s_plot_DOS(obj, coupling)
            if nargin < 2 || isempty(coupling), plot_coupling = false; 
            else plot_coupling = true; end
            
            
            figure
            plot_settings
            
            N_baths = obj.Numel_Baths;
            % how many baths
            if N_baths == 1
                sp1 = 1; sp2 = 1;
            elseif N_baths ==2
                sp1 = 1; sp2 = 2;
            elseif N_baths >2
                sp1 = ceil(N_baths/2); sp2 = 2;
            end
            colors = {[0.5,0.5,1]; [1, 0.5, 0.5]; [0.5, 1, 0.5]};
            for alpha = 1:N_baths
                
                subplot(sp1, sp2, alpha)
                
                % spin yes / no
                if obj.Spin_conserved(alpha)
                    
                    if obj.Spin_symmetric(alpha)
                        Dos =  -imag(diagsum(obj.GF{alpha},1,2)/pi);
                        t = linspace(min(obj.Omega.Points), max(obj.Omega.Points), 500)';
                        Dos_fill = interp1(obj.Omega.Points+obj.Shift_Baths(alpha)*obj.Mu(alpha), Dos,t, 'linear', 'extrap');
                        plot(t,Dos_fill .* obj.Fermi(t,obj.Mu(alpha),obj.Beta(alpha),1))
                        hold on
                        fill([t; max(t); min(t)], [Dos_fill .* obj.Fermi(t,obj.Mu(alpha),obj.Beta(alpha),1); 0;0], colors{alpha})
                        
                        plot(t,Dos_fill)
                        set(gca, 'view', [90,90])
                        ylabel('DOS')
                        xlabel('Energy')
                        
                        
                         spin_text = '$\uparrow = \downarrow$';
                    else
                        Dos_up = imag(diagsum(obj.GF{alpha}(1:obj.Numel_Orbitals(alpha), 1:obj.Numel_Orbitals(alpha),:),1,2))/pi;
                        Dos_down = -imag(diagsum(obj.GF{alpha}(1+obj.Numel_Orbitals(alpha):end, 1+obj.Numel_Orbitals(alpha):end,:),1,2))/pi;
                        t = linspace(min(obj.Omega.Points), max(obj.Omega.Points), 500)';
                        
                        Dos_up_fill = interp1(obj.Omega.Points+obj.Shift_Baths(alpha)*obj.Mu(alpha), Dos_up, t, 'linear', 'extrap');
                        Dos_down_fill = interp1(obj.Omega.Points+obj.Shift_Baths(alpha)*obj.Mu(alpha), Dos_down, t, 'linear', 'extrap');
                        
                        plot(t,Dos_up_fill .* obj.Fermi(t,obj.Mu(alpha),obj.Beta(alpha),1))
                        hold on
                        fill([t; max(t); min(t)], [Dos_up_fill .* obj.Fermi(t,obj.Mu(alpha),obj.Beta(alpha),1); 0;0], colors{alpha})
                        
                        plot(t,Dos_down_fill .* obj.Fermi(t,obj.Mu(alpha),obj.Beta(alpha),1))
                        fill([t; max(t); min(t)], [Dos_down_fill .* obj.Fermi(t,obj.Mu(alpha),obj.Beta(alpha),1); 0;0], colors{alpha})
                        
                        plot(t, Dos_up_fill)
                        set(gca, 'view', [90,90])
                        ylabel('DOS')
                        xlabel('Energy')
                        
                        plot(t, Dos_down_fill)
                        spin_text = '$\uparrow | \downarrow$';
                        ylim([-max(abs(ylim)), max(abs(ylim))]);
                        
                        
                        
                        
                        
                        
                    end
                else
                    Dos = -imag(diagsum(obj.GF{alpha},1,2))/pi;
                    t = linspace(min(obj.Omega.Points), max(obj.Omega.Points), 500)';

                    Dos_fill = interp1(obj.Omega.Points+obj.Shift_Baths(alpha)*obj.Mu(alpha), Dos, t, 'linear', 'extrap');
                    plot(t,Dos_fill .* obj.Fermi(t,obj.Mu(alpha),obj.Beta(alpha),1))
                    
                    hold on
                    fill([t; max(t); min(t)], [Dos_fill .* obj.Fermi(t,obj.Mu(alpha),obj.Beta(alpha),1); 0;0], colors{alpha})
                    plot(t, Dos_fill)
                    set(gca, 'view', [90,90])
                    ylabel('DOS')
                    xlabel('Energy')
                    hold on
                   
                   spin_text = '$\uparrow \downarrow$';
                  
                    
                end
                set(gca, 'xdir', 'reverse')

                set(gca,'xaxislocation', 'origin')
                grid on
                text(max(xlim) - 0.05 * diff(xlim), sum(ylim)/2, spin_text, 'interpreter', 'latex', 'fontsize', 20, 'horizontalalignment', 'center', 'verticalalignment', 'middle');
                if alpha == 2 && obj.Spin_symmetric(alpha)
                    set(gca, 'ydir', 'reverse')
                end
                if plot_coupling
                    
                    
                    %add row and column head lines
                    [C, leng] = print_aligned_table (coupling{alpha});
                    if obj.Spin_symmetric(alpha)
                        spin_labels = repmat('o', 1, obj.Numel_Orbitals(alpha));
                        orbital_numbers = 1:obj.Numel_Orbitals(alpha);
                    else
                        spin_labels = [repmat('u', 1,obj.Numel_Orbitals(alpha)), repmat('d', 1, obj.Numel_Orbitals(alpha))];
                        orbital_numbers = [1:obj.Numel_Orbitals(alpha), 1:obj.Numel_Orbitals(alpha)];
                    end
                    
                    site_numbers = (1:size(coupling{alpha},1))';
                    
                    row = arrayfun(@(x,y) [repmat(' ', 1, floor((leng - length(x)-length(num2str(y)))/2)), x,num2str(y), repmat(' ', 1, ceil((leng - length(x)-length(num2str(y)))/2))], spin_labels, orbital_numbers, 'uniformoutput', false);
                    column = [{repmat(' ', 1, leng)}; arrayfun(@(x) [num2str(x), repmat(' ', 1, leng-2-length(num2str(x))), '| '], site_numbers, 'uniformoutput', false)];

                    D = [column, [row; C]];
                    E = [['Coupling:', repmat(' ', 1, leng*(size(C,2)+1)-9)]; cell2mat(D)];
                    text(max(xlim) - 0.9 * diff(xlim), sum(ylim)/2, E, 'fontsize', 10, 'horizontalalignment', 'center', 'verticalalignment', 'middle');
                end 
                
                %beta and mu
                betamu_text = ['$\beta$: ', num2str(obj.Beta(alpha)), ', $\mu$: ', num2str(obj.Mu(alpha)),', shift: ', num2str(obj.Shift_Baths(alpha))];
                text(max(xlim) - 0.1 * diff(xlim), sum(ylim)/2, betamu_text, 'interpreter', 'latex', 'fontsize', 15, 'horizontalalignment', 'center', 'verticalalignment', 'middle');
               
                
            end
            
        end
        
        
        
       
        
        function values = f_G_eval(obj, x, alpha, c_spin, row, column)
            if nargin < 6 || isempty(column), column = 1; end
            if nargin < 5 || isempty(row), row = 1; end
            if nargin < 4 || isempty(c_spin), c_spin = -1; end
            if nargin < 3 || isempty(alpha), alpha = 1; end
            
            if numel(obj.GF{alpha,1}) == 1
                values = obj.GF{alpha,1};
            else
            
                pp_G = spline(obj.Omega.Points, squeeze(obj.GF{alpha,1}(row,column,:)));
            
                values = ppval(pp_G, x); 
            end
        end
        
        function num_orb = Numel_Orbital(obj, bath_id, spin_id)
            %NOTE: help function to use new definition of Numel_Orbital (for transfere to class "Bath"
            error('use Numel_Orbitals now')
            if nargin < 3 || isempty(spin_id), spin_id = 1; end
            if nargin < 2 || isempty(bath_id), bath_id = 1:obj.Numel_Baths; end
            num_orb = sqrt(obj.Numel_Orbitals(bath_id, spin_id));
        end
            
            
    end
    
    methods(Static = true) 
        
        
         function f_check_Coupling(bath, Coupling, basis)
             if isa(bath, 'Bath')
                 numel_orbitals = bath.Numel_Orbitals; 
                 numel_baths = bath.Numel_Baths;
                 spin_symmetric = bath.Spin_symmetric;
                 spin_conserved = bath.Spin_conserved;
                 GF = bath.GF;
                 n_leads = bath.Numel_Baths;
             else
                 numel_orbitals = bath.numel_orbitals; 
                 numel_baths = bath.N_leads;
                 spin_symmetric = bath.spin_symmetric; 
                 spin_conserved = bath.spin_conserved;
                 GF = bath.GF;
                 n_leads = bath.N_leads;
             end
            
             
             if ~isa(bath,'Bath') %NOTE t_leads is just a property of bath in setup if tight binding is chosen
                 if strcmp(bath.type, 'from_file')
                     for k = 1:n_leads
                         if size(GF{k},1) ~= numel_orbitals(k)
                             error(['Numel of orbitals (', num2str(numel_orbitals(k)), ') does not match dimension of input Green function (', num2str(size(GF{k},1)), ')'])
                         end
                     end
                 elseif strcmp(bath.type, 'tight_binding')
                     for k = 1:n_leads
                         if size(bath.t_leads{k},1) ~= numel_orbitals(k)
                             error(['Numel of orbitals (', num2str(numel_orbitals(k)), ') does not match dimension of input t leads (', num2str(size(bath.t_leads{k},1)), ')'])
                         end
                     end
                 end
             else
                 for k = 1:n_leads
                     if size(GF{k},1) ~= numel_orbitals(k)
                         error(['Numel of orbitals (', num2str(numel_orbitals(k)), ') does not match dimension of Green function (', num2str(size(GF{k},1)), ')'])
                     end
                 end
             end
             
            % ---- check dimensions of input variables
            [siz_Coupl1, siz_Coupl2] = size(Coupling); %n_bath, n_spin
            n_bath_orbitals = cellfun(@(x) size(x,2), Coupling);
            n_system_sites = cellfun(@(x) size(x,1), Coupling);
            
            %{num_baths, num_spin}(n_orbitals x n_orbitals)
            
            if siz_Coupl1 ~= numel_baths
                error(['Number of Baths in Coupling variable does not match bath (#cell rows), n_baths: ', num2str(obj.Numel_Baths)]);
            end
%             if siz_Coupl2 ~= size(obj.GF,2) %TODO: maybe define as constant in Bath, when initialized
%                 error(['Number of Spins in Coupling variable does not match bath (#cell columns), n_spins',  num2str(size(obj.GF,2))])
%             end

            %Note: the size of the Green functions is twice the size of the number of orbitals if there are
            %differences in the spin!
            if any(n_bath_orbitals(:)./(2-spin_symmetric) ~= (numel_orbitals(:))) 
                error(['Number of Orbitals in Coupling variable does not match bath, #bath orbitals: ', num2str(numel_orbitals')])
            end
            
            if any(basis.N_site ~= n_system_sites)
                error(['Number of system sites in Coupling variable does not match basis of central system: ', num2str(basis.N_site)])
            end
            
            disp('---------------coupling check okay :) -------------------')
            
       
            
            
        end
        
        function G_TB = sf_tight_binding_single_GF(grid, t_leads, eps_0)
            if nargin < 3 || isempty(eps_0), eps_0 = 0; end
           
                %cell array of Green functions (bath and spin)
            %numel(orb) size(t_leads)
            G = @(om,eps_0, t) (om-eps_0)/(2*t^2) - sign(om-eps_0 + 2*t).*sqrt((om - eps_0).^2/(4*t^4) - 1/t^2);
            % returns a Green function object retained by approximation of 
            G_TB = G(grid.Points, eps_0, t_leads);
        end
        
        function G = sf_dos_2_GR(grid, dos)
            %TODO: Think about integration boudaries! - put into integration function
            warning('Use sf_GR_imag_2_real(grid, -pi*1i*real(dos)')
            %NOTE: if grid is a class or not; think of general way!
            G = Bath.sf_GR_imag_2_real(grid, -pi*1i*real(dos));
            
%             if ~isa(grid, 'Grid')
%                 omega = grid;
%                 i0plus = 10^-4;
%             else
%                 omega = grid.Points;
%                 i0plus = grid.Imag_Infinitesimal;
%             end
%                         
%             G_imag = real(dos) *-pi*1i;
%             G_pp = spline(omega, G_imag);
%             G_real = -1/pi *real(integral(@(x) imag(ppval(G_pp,x))./(omega-x+i0plus*1i),-10,10,'arrayvalued', true));
%             
%             G = zeros(1,1,numel(omega));
%             %TODO: think of Green function class, Brainstorm which features are necessary
%             G(1:end) = G_real + 1i*G_imag;
        end
            
        function G_real = sf_GR_imag_2_real(grid, G_imag)
            omega = grid.Points;
            i0plus = grid.Imag_Infinitesimal;
            
            %NOTE: works also efficiently with blockdiagonal Green functions since 0 GF are ignored
            
            %NOTE: effective trapezoidal integration - normal trapz for vectors end up with memory overflow... -> workaround
            % use of v * v' is matrix (new in Matlab)
            G_real = zeros(size(G_imag));
            for j = 1:size(G_imag,1)
                for k = 1:size(G_imag,2)
                    if j == 1 && k == 1
                        tic
                    end
                    GI = squeeze(G_imag(j,k,:));
                    
                    if sum(abs(GI)) < 10^-10
                        G_real(j,k,:) = zeros(numel(GI),1);
                    else
                        
                        
                        if grid.N > 5000
                            intervals = 0:5000:grid.N;
                            if intervals(end) < grid.N, intervals = [intervals, grid.N]; end
                            for r = 1:numel(intervals)-1
                                interval = (1+intervals(r)): intervals(r+1);
                                G_real(j,k,interval) = 1/pi*real(sum((omega((2:end)) - omega((1:end-1))).*...
                                    (imag(GI((2:end)))./((omega(interval(1:end)))'-omega((2:end)) + i0plus*1i) + imag(GI((1:end-1)))./((omega(interval(1:end))')-omega((1:end-1))+ i0plus*1i)),1)/2).';
                            end
                        else
                            
                            
                            G_real(j,k,:) = 1/pi*real(sum((omega(2:end) - omega(1:end-1)).*(imag(GI(2:end))./((omega(1:end))'-omega(2:end) + i0plus*1i) + imag(GI(1:end-1))./((omega(1:end)')-omega(1:end-1)+ i0plus*1i)),1)/2).';
                        end
                    end
                    if j == 1 && k == 1
                        time = toc;
                        disp(['Expected time for calculation of Greens function: ', num2str(time), ' s'])
                    end
                end
            end
            %plot(omega, G_real)
            
            %NOTE: vector valued integral quite expansive - compare with trapz! poles could be treated extra!
            
            if false
                G_pp = spline(omega, G_imag);
                G_real = -1/pi *real(integral(@(x) imag(ppval(G_pp,x))./(omega-x+i0plus*1i),grid.Limits(1),grid.Limits(2),'arrayvalued', true));
            end
        end
        
        function G_imag = sf_GR_real_2_imag(grid, G_real)
            omega = grid.Points;
            i0plus = grid.Imag_Infinitesimal;
            %TODO: convergence may not be given, if limits of integral are too limited, think about version where
            %it is possible to estimate those!
            
            G_imag = zeros(size(G_real));
            %loop over entries
            for j = 1:size(G_real,1)
                for k = 1:size(G_real,2)
                    if j == 1 && k == 1
                        tic
                    end
                    GR = squeeze(G_real(j,k,:));
                    
                    if sum(abs(GR)) < 10^-10
                        G_real(j,k,:) = zeros(numel(GI),1);
                    else
                        if grid.N > 5000
                            intervals = 0:5000:grid.N;
                            if intervals(end) < grid.N, intervals = [intervals, grid.N]; end
                            for r = 1:numel(intervals)-1
                                interval = (1+intervals(r)): intervals(r+1);
                                G_imag(j,k,interval) = -1/pi*real(sum((omega((2:end)) - omega((1:end-1))).*...
                                   (real(GR((2:end)))./((omega(interval(1:end)))'-omega((2:end)) + i0plus*1i) + real(GR((1:end-1)))./((omega(interval(1:end))')-omega((1:end-1))+ i0plus*1i)),1)/2).';
                                
                               
                            end
                        else
                        
                            %G_pp = spline(omega, squeeze(G_real(j,k,:)));
                            G_imag(j,k,:) = 1/pi * real(sum((omega(2:end) - omega(1:end-1)).*(real(GR(2:end))./((omega(1:end))'-omega(2:end) + i0plus*1i) + real(GR(1:end-1))./((omega(1:end)')-omega(1:end-1)+ i0plus*1i)),1)/2).';
                        end
                    end
                    if j == 1 && k == 1
                        time = toc;
                        disp(['Expected time for calculation of Greens function: ', num2str(time), ' s'])
                    end
                end
            end
            
            
                    
            %G_pp = spline(omega, G_real);
            %G_imag =  -1/pi *real(integral(@(x) real(ppval(G_pp,x))./(omega-x+i0plus*1i),-10,10,'arrayvalued', true));
            if (G_imag(1) > 10^-8) || (G_imag(end) > 10^-8), warning('Convergence Problems - Integration area may not be sufficient'); end
            
        end
        
        function G_flatband = sf_flatband_from_GR(grid, G, eps_0, reconstruct_real_part)
            %flatband estimation given a retarded Green function (rather Tight Binding)
            %find G at eps_0!
            if nargin < 3 || isempty(eps_0), eps_0 = 0; end
            if nargin < 4 || isempty(reconstruct_real_part), reconstruct_real_part = false; end
            
            G_flatband_height = imag(Bath.sf_evaluate_G(grid, G, eps_0));
            
            G_flatband_imag = permute((abs(grid.Points- eps_0 ) < pi/(2*abs(G_flatband_height))) * 1i*G_flatband_height, [3,2,1]);
            if reconstruct_real_part
                G_flatband_real = Bath.sf_GR_imag_2_real(grid, G_flatband_imag);
            else
                G_flatband_real = 0;
            end
            G_flatband = G_flatband_real + G_flatband_imag;
            
            
        end
        
        function values =  sf_evaluate_G(grid, G, x)
            %evaluate Green function at point x
            pp_G = spline(grid.Points, G);
            
            values = ppval(pp_G, x);
        end
            
            
        function bath = c_tight_binding(GF_grid, t_leads, shift_bath, eps_0, integrator, N_spin, N_orbitals, beta, mu)
            if nargin < 8, beta = []; end
            if nargin < 9, mu = []; end
            if nargin < 5 || isempty(integrator), integrator = 'Fast'; end
            if nargin < 4 || isempty(eps_0), eps_0 = zeros(size(t_leads)); end
            if nargin < 3 || isempty(shift_bath), shift_bath = false(size(t_leads)); end
            
            
%             %test: 
%             t_leads = {[2 0.1; 0.1, 3]; 3};
%             shift_bath = [true; false];
%             eps_0 = {[-1, 0; 0, 0.3]; 0};
%             
%             integrator = 'Fast';
%             N_spin = 2; 
%             N_orbitals = [2;1];
%             
%             %test2:
%             t_leads = {[2 0; 0 1]; [3 0; 0 1]};
%             shift_bath = [true; true];
%             eps_0 = {[-1, 0 ; 0 1]; [-1, 0 ; 0 1]};
%             N_orbitals = [1;1];
%             integrator = 'FastDetail';
%             N_spin = 2;
            
            
            %GR{bath}(orbital_i, orbital_j, omega)
            if iscell(t_leads)  %NOTE: new procedure for many orbitals in one bath / lead
                %NOTE: cell array indicates the number of baths (a different spin could be an indication for a
                %different bath
                %the number of baths is indicated by the first dimension of the columns, the second dimension
                %indicates a difference in spin (spin could also be implemented in the orbitals
                
                N_baths = size(t_leads,1);
                Spin_symmetric = false(N_baths,1);
                Spin_conserved = false(N_baths,1);
                if N_spin == 1
                    % spinless
                   
                    Spin_symmetric(:,1) = true(N_baths,1);
                    Spin_conserved(:,1) = true(N_baths,1);
                    GF_ret = cell(N_baths, 1);
                    for alpha = 1:N_baths
                        N_orb = N_orbitals(alpha);
                        G_temp = zeros(N_orb, N_orb, GF_grid.N);
                        
                        V = tril(t_leads{alpha,1});
                        s = zeros(size(t_leads{alpha,1}));
                        tL = eps_0{alpha,1} + (t_leads{alpha,1}-diag(diag(t_leads{alpha,1})));
                        S = 1.*eye(size(t_leads{alpha,1}));
                        omega_range = GF_grid.Points;
                        %Sancho parameters:
                        NullPlus = 4;
                        epsilon = 10^-7;
                        nmax = 1000;
                        
                        for omega_ind = 1:numel(omega_range)
                            omega = omega_range(omega_ind);
                            G_temp(:,:,omega_ind)  = Bath.sf_semi_chain_sancho( V,s,tL,S,omega,NullPlus,epsilon,nmax );
                        end
                        
                        GF_ret{alpha,1} = G_temp;
                    end
                    
                else
                    
                    
                   
                    if size(eps_0,1) ~= N_baths || size(shift_bath,1) ~= N_baths || size(N_orbitals,1) ~= N_baths
                        error('in construction of c_tight_binding bath: t_leads has to match eps_0')
                    end
                    
                    
                    
                    
                    GF_ret = cell(N_baths, 1);
                    for alpha = 1:N_baths
                        N_orb = N_orbitals(alpha);
                        [nx] = size(t_leads{alpha},1);
                        if nx == N_orbitals(alpha,1)
                            % spinsymmetric and spin conserved
                            Spin_symmetric(alpha,1) = true;
                            Spin_conserved(alpha,1) = true;
                            G_temp = zeros(N_orb, N_orb, GF_grid.N);
                            
                            
                            V = tril(t_leads{alpha,1});
                            s = zeros(size(t_leads{alpha,1}));
                            tL = eps_0{alpha,1} + (t_leads{alpha,1}-diag(diag(t_leads{alpha,1})));
                            S = 1.*eye(size(t_leads{alpha,1}));
                            omega_range = GF_grid.Points;
                            %Sancho parameters:
                            NullPlus = 4;
                            epsilon = 10^-7;
                            nmax = 1000;
                            
                            for omega_ind = 1:numel(omega_range)
                                omega = omega_range(omega_ind);
                                G_temp(:,:,omega_ind)  = Bath.sf_semi_chain_sancho( V,s,tL,S,omega,NullPlus,epsilon,nmax );
                            end
                            
                            
                        elseif nx == 2*N_orbitals(alpha,1)
                            Spin_symmetric(alpha,1) = false;                                
                            G_temp = zeros(2*N_orb, 2*N_orb, GF_grid.N);

                            if all(all(t_leads{alpha}(1:N_orb,1+N_orb:2*N_orb) == 0)) && all(all(t_leads{alpha}(1+N_orb:2*N_orb, 1:N_orb) == 0))
                                Spin_conserved(alpha,1) = true;

                                
                                %spin up
                                spin_up = 1:N_orb;
                                V = tril(t_leads{alpha,1}(spin_up, spin_up));
                                s = zeros(size(t_leads{alpha,1}(spin_up, spin_up)));
                                tL = eps_0{alpha,1}(spin_up, spin_up) + (t_leads{alpha,1}(spin_up, spin_up)-diag(diag(t_leads{alpha,1}(spin_up, spin_up))));
                                S = 1.*eye(size(t_leads{alpha,1}(spin_up, spin_up)));
                                omega_range = GF_grid.Points;
                                %Sancho parameters:
                                NullPlus = 4;
                                epsilon = 10^-7;
                                nmax = 1000;
                                
                                if all(V(:) == 0)
                                    G_temp(spin_up, spin_up, :) = zeros(numel(spin_up), numel(spin_up), numel(omega_range));
                                else
                                    
                                    for omega_ind = 1:numel(omega_range)
                                        omega = omega_range(omega_ind);
                                        G_temp(spin_up, spin_up,omega_ind)  = Bath.sf_semi_chain_sancho( V,s,tL,S,omega,NullPlus,epsilon,nmax );
                                    end
                                end
                                
                                %spin down
                                spin_down = 1+N_orb:2*N_orb;
                                V = tril(t_leads{alpha,1}(spin_down, spin_down));
                                s = zeros(size(t_leads{alpha,1}(spin_down, spin_down)));
                                tL = eps_0{alpha,1}(spin_down, spin_down) + (t_leads{alpha,1}(spin_down, spin_down)-diag(diag(t_leads{alpha,1}(spin_down, spin_down))));
                                S = 1.*eye(size(t_leads{alpha,1}(spin_down, spin_down)));
                                omega_range = GF_grid.Points;
                                %Sancho parameters:
                                NullPlus = 4;
                                epsilon = 10^-7;
                                nmax = 1000;
                                if all(V(:) == 0)
                                    G_temp(spin_down, spin_down, :) = zeros(numel(spin_down), numel(spin_down), numel(omega_range));
                                else
                                    for omega_ind = 1:numel(omega_range)
                                        omega = omega_range(omega_ind);
                                        G_temp(spin_down, spin_down,omega_ind)  = Bath.sf_semi_chain_sancho( V,s,tL,S,omega,NullPlus,epsilon,nmax );
                                    end
                                end
                                
                                
                            else
                                Spin_conserved(alpha,1) = false;
                                
                                V = diag(diag(t_leads{alpha,1}));
                                % V = tril(t_leads{alpha,1}); %crosshoppings
                                % O -a- O --- O
                                % | c / |     |
                                % b  X  |
                                % | d \ |
                                % O --- O --- O
                                % V = [a,c; d a]; tL = [eps_1, b; b, eps_2]
                                s = zeros(size(t_leads{alpha,1}));
                                tL = eps_0{alpha,1} + (t_leads{alpha,1}-diag(diag(t_leads{alpha,1})));
                                S = 1.*eye(size(t_leads{alpha,1}));
                                omega_range = GF_grid.Points;
                                %Sancho parameters:
                                NullPlus = 4;
                                epsilon = 10^-7;
                                nmax = 1000;
                                
                                for omega_ind = 1:numel(omega_range)
                                    omega = omega_range(omega_ind);
                                    G_temp(:, :,omega_ind)  = Bath.sf_semi_chain_sancho( V,s,tL,S,omega,NullPlus,epsilon,nmax );
                                end   
                               
                            end
                        end
                            
                        GF_ret{alpha,1} = G_temp;
                    end
                    
                end
                
            else
                error('input must be given as cell - new bath class')
                
                %left lead is negative!
            end
            bath = Bath(GF_ret, GF_grid, shift_bath, integrator, N_spin, N_orbitals, Spin_conserved, Spin_symmetric, beta, mu);
            
        end
        
        
        function bath = c_constant(grid, N_leads)
            gf = num2cell(-1i*1*ones(N_leads, 1));
            integrator = 'Constant';
            shift_bath = false;
            
            bath = Bath(gf, grid, shift_bath, integrator);
        end
      
        
        function bath = c_wide_band(GF_grid, t_leads, shift_bath, eps_0, cutoff, reconstruct_real_part, integrator)
            %NOTE: no orbitals! no spin difference (spin symmetric and spin conserved == spinless)
            %details:
            % for each lead there can be a unique t_leads, and shift and eps_0, cutoff can also be uniquly defined
            % t_leads ... defines the height and the width (in case of cutoff = true) of the wide band für each bath
            % shift   ... defines if the lead is shifted in combination with mu_\alpha (just relevant for cutoff = true case
            % eps_0   ... center of bandwith (if cutoff = true)
            % cutoff  ... defines a cutoff, will create a Green's function, cutoff = false creates a single entry
            % in Green's function -> evaluation is then faster
            if nargin < 7 || isempty(integrator), integrator = 'Fast'; end
            if nargin < 6 || isempty(reconstruct_real_part), reconstruct_real_part = true; end
            if nargin < 5 || isempty(cutoff), cutoff = true; end
            N_baths = size(t_leads,1);
            
            if nargin < 4 || isempty(eps_0), eps_0 = cell(N_baths,1); [eps_0{:}] = deal(0); end
            if nargin < 3 || isempty(shift_bath), shift_bath = false; end
            
            %test 
%             t_leads = {2;3};
%             eps_0 = {-1;0};
%             
            
            
            GF = cell(N_baths,1);
            N_orbitals = ones(N_baths,1);
            Spin_conserved = true(N_baths,1);
            Spin_symmetric = true(N_baths,1);
            N_spin = 2; %makes no difference if spinless and spin_conserved + spin_symmetric
            G_flatband_height = zeros(N_baths,1);
            

            
            
            for alpha = 1:N_baths
                G_flatband_height(alpha) = 1./abs(t_leads{alpha});
                %obj.Cutoff = eps_0 +t_leads*pi/2 *[-1,1];
                
                
                if cutoff == true
                    
                    
                    G_flatband_imag = permute((abs(GF_grid.Points - eps_0{alpha,1}) < pi/2*abs(t_leads{alpha,1}) )* -1i*G_flatband_height(alpha), [3,2,1]);
                    
                    if reconstruct_real_part
                        G_flatband_real = Bath.sf_GR_imag_2_real(GF_grid, G_flatband_imag);
                    else
                        G_flatband_real = 0;
                    end
                    GF{alpha} = G_flatband_real + G_flatband_imag;
                    
                    
                else
                    GF = num2cell(-1i*G_flatband_height);
                end
                
            end
            
           
            
            bath = Bath(GF, GF_grid, shift_bath, integrator, N_spin, N_orbitals, Spin_conserved, Spin_symmetric);

            
        end
        
        function [ G_semi ] = sf_semi_chain_sancho( V,s,tL,S,omega,NullPlus,epsilon,nmax )
            
            %
            %   Sancho-Methode nach "Journal of Physics F: Metal Physics 14 (5), 1205
            %   (1984)".
            %
            % by Michael Rumetshofer
            if all(V(:) == 0)
                warning('Sancho does not work well in the atomic limit (V = 0)')
            end
            % Initial values.
            epss  = S*((omega)+1i*10^(-NullPlus)) - tL;   %omega + i0+ -tL
            eps   = S*((omega)+1i*10^(-NullPlus)) - tL;
            alpha = -(s'*((omega)+1i*10^(-NullPlus)) - V');  % = V' (if s == 0)
            beta  = -(s*((omega)+1i*10^(-NullPlus)) - V);    % = V  (if s == 0)
            
            % Convergence.
            converge  = 0;
            n         = 0;
            
            %
            while(~converge)
                a = eps\alpha;
                b = eps\beta;
                epss = epss - alpha*b;
                eps  = eps - beta*a - alpha*b;
                
                % Proof if convergence of T fulfilled.
                if (sum(sum(((abs(alpha*beta))<= epsilon))) == numel(alpha))
                    converge = 1;
                end
                if (n>nmax)
                    converge = 1;
                    fprintf('\nG semi chain is not converged!')
                end
                
                alpha = alpha*a;
                beta  = beta*b;
                
                n      = n + 1;
            end
            
            G_semi = inv(epss);
            
        end
        
        
        
    end
        
             
        
    
    
end

