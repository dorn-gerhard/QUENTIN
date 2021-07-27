classdef Bath
    %author: Gerhard Dorn
    %date: December 2018
    %use: describes the (non-interacting) bath including orbitals used in the master equation approach, bath is described by Temperature (beta), Energy
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
    
    
    
    % Define what Bath must do (in the most general case)
    % Then distinguish easier cases and implement simplifications!
    
    properties
        N_spin
        Spin_conserved
        Spin_symmetric
        Numel_Orbitals
        GF
        
        Integrator
        Shift_Bath
        
        
        %use a precalc object???
        
        
    end
    
    
    methods
        function obj = Bath(GF, Integrator, Shift_Bath, N_spin, Spin_conserved, Spin_symmetric, Numel_Orbitals)
            
            if isa(GF, 'Green')
                obj.GF = GF;
            else
                error('Input for Bath class shall be a Green function (class Green)')
            end
            
            obj.Integrator = Integrator;
            obj.Shift_Bath = Shift_Bath;
            obj.N_spin = N_spin;
            obj.Spin_conserved = Spin_conserved;
            obj.Spin_symmetric = Spin_symmetric;
            obj.Numel_Orbitals = Numel_Orbitals;
            
            
            %check if dimensions fit
            if obj.GF.Dim1 ~= obj.Numel_Orbitals
            
            end
        
        
        end
        
        
        function GR = f_GR_shifted(obj, shift)
            %returns the GF shifted by the value mu: G(omega - shift)
            
          
            if ~obj.Shift_Bath || shift == 0 
                GR = obj.GF.R;
            else
                [nx,ny,nz] = size(obj.GF.R);
                GR_extended = zeros(nx,ny,nz+1);
                if shift > 0
                    grid_extended = [obj.GF.Grid.Limits(1)- shift; obj.GF.Grid.Points];
                    GR_extended(1:nx,1:ny,2:end) = obj.GF.R;
                    y1 = obj.GF.R(1:nx,1:ny,1);
                    y2 = obj.GF.R(1:nx,1:ny,2);
                    x1 = obj.GF.Grid.Points(1);
                    x2 = obj.GF.Grid.Points(2);
                    %linear extrapolation
                    GR_extended(1:nx,1:ny,1) = ( (y2-y1) * grid_extended(1) + y1*x2-y2*x1)/(x2-x1);
                elseif shift < 0
                    grid_extended = [obj.GF.Grid.Points; obj.GF.Grid.Limits(2) - shift]; %mu is negative
                    GR_extended(1:nx,1:ny,1:end-1) = obj.GF.R;
                    %linear extrapolation
                    y1 = obj.GF.R(1:nx,1:ny,end-1);
                    y2 = obj.GF.R(1:nx,1:ny,end);
                    x1 = obj.GF.Grid.Points(end-1);
                    x2 = obj.GF.Grid.Points(end);
                    GR_extended(1:nx,1:ny,end) = ( (y2-y1) * grid_extended(end) + y1*x2-y2*x1)/(x2-x1);
                end
                
                %GR = reshape(spline(grid_extended, GR_extended, obj.GF.Grid.Points - shift, 'linear'),nx,ny,nz);
                %linear interpolation better suited
                GR = reshape(interp1(grid_extended, reshape(GR_extended, nx*ny,[]).', obj.GF.Grid.Points - shift, 'linear').',nx,ny,nz);
                %plot(obj.GF.Grid.Points, squeeze(real(GR(1,1,:))))
                %hold on
                %plot(obj.GF.Grid.Points, squeeze(real(obj.GF.R(1,1,:))))
            end
        end
        
        function [C,obj,warn_precalc] = f_lookup(obj)
            %TODO: rethink and make property of equilibrium class, 
            %respectively: rethink for Non-equilibrium classes!!!
        end
        
        
        function F_calc = f_calculate(obj)
            % shall work for all situations, start with non-equilibrium, simplify for equilibrium
            % implement different integration techniques
            % implement storage and interpolation
            
        end
        
        
        function Spin_id = Spin_fun(obj,x)
            warning('check if this function is really needed')
            if obj.N_Spin == 1
                Spin_id =  1;
            elseif obj.N_Spin == 2
                Spin_id = x/2+1.5;
            else
                error('unusual number of spins! Not implemented!')
            end
        end
    end
    
    methods (Static =true)
        
            
    end
    
    
end