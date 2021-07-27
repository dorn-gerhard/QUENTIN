classdef Green
    
    %author: Gerhard Dorn
    %date: january 17
    %use: handle Green functions
    %features: 
    %          handle:
    %             - basic operations such as dagger, multiplication, inversion
    %             - storage of retarded and keldysh part
    %          analyse:
    %
    %          initialize:
    
    properties
        R           % retarded part
        K           % keldysh part
        Grid        % Grid
        Dim1        % first dimension
        Dim2        % second dimension
    end
    
    
    methods
        function obj = Green(grid, GR, GK, s_func)
            % constructor
            
            %NOTE: ensure that the paths of the two functions mmat and multinv are added
            
            if nargin < 3 || isempty(GK), GK = zeros(size(GR)); end
            if nargin < 1 || isempty(grid), error('Definition of Grid necessary!'); end
            if nargin >= 4 && ~isempty(s_func) 
                if ~isdeployed()
                    disp('Keldysh part calculated using s_func')
                end
                GK = Green.f_keldysh_from_equilibrium(GR, s_func);
            end
         
            
            obj.R = GR;
            obj.K = GK;
            obj.Grid = grid;
            obj.Dim1 = size(GR,1);
            obj.Dim2 = size(GR,2);
            
            %checks: 
            if  ~all(size(GR) == size(GK)) 
                error('Sizes of G_ret and G_K do not match!');
            end
            if size(GR,3) ~= grid.N
                error('G_ret does not match Grid points!')
            end
            
        end
        
        function GA = A(obj)
            
            if obj.Dim1 ~= obj.Dim2
                GA = obj.R;
            else
                GA = Green.f_dagger(obj.R);
            end
           
        end
        

        
        function GK = f_keldysh_inverse(obj)
            % -R^-1 * K * A^-1
            if obj.Dim1 ~= obj.Dim2, error('inverse only for square Green functions defined!'); end
            
            switch obj.Dim1
                case 1
                    
                    temp_K = -1./obj.R .* obj.K ./ conj(obj.R);
                    
                otherwise
                    R_inv = Green.f_inverse(obj.R);
                    
                    A_inv = Green.f_dagger(R_inv);
                    
                    temp_K = -mmat(mmat(R_inv, obj.K), A_inv);
            
                    
            end
            
            GK = temp_K;
            
        end
        
        
        function obj = f_copy_spin(obj)
            
            obj.K = [obj.K, zeros(size(obj.K)); zeros(size(obj.K)), obj.K];
            obj.R = [obj.R, zeros(size(obj.R)); zeros(size(obj.R)), obj.R];
            obj.Dim1 = obj.Dim1*2;
            obj.Dim2 = obj.Dim2*2;
        end
        
        function s_plot(obj, c_site_1, c_site_2)
           
            
            if nargin < 2 || isempty(c_site_1), c_site_1 = 1; end
            if nargin < 3 || isempty(c_site_2), c_site_2 = 1; end
             figure
            k = c_site_1;
                
            j = c_site_2;
            
            if ~isempty(obj.K)
                use_k = 2; else 
                use_k = 1; 
            end
            
            
            
            plot_settings
           
            subplot(2, use_k,1)
            plot(obj.Grid.Points, real(squeeze(obj.R(k,j, :))))
            title('Retarded, real')
            subplot(2, use_k, 1 + use_k)
            plot(obj.Grid.Points, imag(squeeze(obj.R(k,j, :))))
            title('Retarded, imag')
            
            if use_k == 2
                subplot(2, use_k,2)
                plot(obj.Grid.Points, real(squeeze(obj.K(k,j, :))))
                title('Keldysh, real')
                subplot(2, use_k, 4)
                plot(obj.Grid.Points, imag(squeeze(obj.K(k,j, :))))
                title('Keldysh, imag')
            end
            plot_supertitle(['$G_{',num2str(c_site_1), num2str(c_site_2), '}$, dim1: ', num2str(obj.Dim1), ', dim2: ', num2str(obj.Dim2)])
            
        end
        
        function total_spectral = s_plot_spectral(obj)
            
            figure
            plot_settings
            
            subplot(2,1,1)
            legend_text = cell(1,obj.Dim1);
            total_spectral = zeros(obj.Grid.N,1);
            for k =1:obj.Dim1
                
                
                
                
                plot(obj.Grid.Points, -1/pi*imag(squeeze(obj.R(k,k, :))))
                hold on
                legend_text{k} =   ['$A_{',num2str(k), num2str(k), '}$'];
                total_spectral = total_spectral + -1/pi*imag(squeeze(obj.R(k,k,:)));
            end
            legend(legend_text)
            
            subplot(2,1,2)
            title('Total spectral function')
            plot(obj.Grid.Points, total_spectral)
            plot_supertitle(['Spectral function, dim1: ', num2str(obj.Dim1), ', dim2: ', num2str(obj.Dim2)])
            
            
        end
        

        
        
    end
    
    methods (Static)
        
        function GF = f_times(A,B)
            if isa(A,'Green') || isa(B, 'Green')
                error('f_times need matrices (only retarded part) as input');
            end
            GF = mmat(A,B);
        end
        
        function A = f_add(A,B)
            if (A.Dim1 ~= B.Dim1) || (A.Dim2 ~= B.Dim2)
                error('Dimensions must match!'); 
            end
            A.R = A.R + B.R;
            A.K = A.K + B.K;
        end
        
        function A = f_minus(A,B)
            if (A.Dim1 ~= B.Dim1) || (A.Dim2 ~= B.Dim2)
                error('Dimensions must match!'); 
            end
            A.R = A.R - B.R;
            A.K = A.K - B.K;
        end
        
        function GF = f_full_inverse(GF)
            
            %retarded
            
            
            temp_R = Green.f_inverse(GF.R);
            temp_K = GF.f_keldysh_inverse();
            
            GF.R = temp_R;
            GF.K = temp_K;
            
            
        end
        
        function GF = f_flip(GF)
            GF.R = permute(GF.R,[2,1,3]);
            GF.K = permute(GF.K,[2,1,3]);
            GF.Dim1 = size(GF.R,1);
            GF.Dim2 = size(GF.R,2);
        end
        
        
        function GF = f_inverse(A)
            
            dim = size(A);
            if dim(1) ~= dim(2)
                error('inverse only for square Green functions defined!')
            end
            
            switch dim(1)
                case 1
                    GF = 1./A;
                case 2
                    L_det =   A(1,1,:) .* A(2,2,:) - A(1,2,:) .* A(2,1,:);
                    GF(1,1,:) = A(2,2,:)./L_det;
                    GF(2,2,:) = A(1,1,:)./L_det;
                    GF(1,2,:) = A(1,2,:)./L_det;
                    GF(2,1,:) = A(2,1,:)./L_det;
                otherwise
                    GF = multinv(A);
            
            end
            
            
        end
        
        
        
        
        function GF = f_dagger(A)
            
            GF = conj(permute(A,[2,1,3]));
        end
        
        function GK = f_keldysh_from_equilibrium(GR, s_func)
            if ndims(s_func) == 2
                s_func = reshape(s_func,1,1,[]); 
            end
            if numel(s_func) ~= size(GR,3)
                error('s_func not adopted to grid!')
            end
            %NOTE: for equilibrium situation s_func(omega) =  1 - 2* fermi_function(omega, mu, beta),
            %NOTE: s_func can be altered / adopted to an effective s_func
            
            
            GK = repmat(s_func, size(GR,1), size(GR,2), 1) .* (GR - Green.f_dagger(GR));
        end
        
        function GK = f_times_full(x,y)
            if ~isa(x,'Green') || ~isa(y,'Green')
                error('both inputs must be Green classes, use f_times for matrices')
            end
            %NOTE: A and B should be Green functions
            if x.Dim2 ~= y.Dim1
                error('dimensions for f_times multiplication must match!')
            end
            R_temp = Green.f_times(x.R, y.R);
            
            K_temp = Green.f_times(x.R, y.K) + Green.f_times(x.K, y.A);
            
            %NOTE: recycle x to return
            
            x.R = R_temp;
            x.K = K_temp;
            
            %correct dimension:
            x.Dim2 = y.Dim2;
            GK = x;
        end
        
      
            
        function [Ir, Ik] = f_integrate_full(green_function, points)
            %integrates over the input, which can either be a Green class or a 3-dim Green function
            
            %TODO: define exponential decay at the end of the grid (via fitting) and add to the integration
            
            
            if isa(green_function, 'Green')
                if nargin == 2
                    warning('The integration points are taken from the Green class')
                end
                points = green_function.Grid.Points;
                Ir = permute(trapz(points, permute(green_function.R,[3,1,2])), [2,3,1]);
                Ik = permute(trapz(points, permute(green_function.K,[3,1,2])), [2,3,1]);
            else
            
                Ir = permute(trapz(points, permute(green_function,[3,1,2])), [2,3,1]);
                Ik = [];
            end
            
            
            % idea of integration of the rest
            
            
        end
        
        function G_full = f_dyson(G_local, G_self_energy)
            
            
            
            E = repmat(eye(G_local.Dim1, G_local.Dim2), 1,1, G_local.Grid.N);
            
            
            Gr = Green.f_times(Green.f_inverse(E - Green.f_times(G_local.R, G_self_energy.R)), G_local.R);
            
            U = Green.f_inverse(E - Green.f_times(G_local.R, G_self_energy.R));
            
            
            Gk = -Green.f_times(Green.f_times(U, ...
                (G_local.K - Green.f_times(G_local.R, Green.f_times(G_self_energy.K, G_local.A)))), ...
                Green.f_dagger(U));
            G_full = Green(G_local.Grid, Gr, Gk);
        end

        
        

        
        
    end
    
    
end
            