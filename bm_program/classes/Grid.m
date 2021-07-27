classdef Grid
    
    %author: Gerhard Dorn
    %date: September 16
    %use: handles Grid
    % 
    %features: 
    %          handle:
    %             - stores a grid, not necessarily uniform
    %             - conversion to other grids
    %             - interpolation between grid points
    %             
    %          visualize:
    %             - print / visualize bath / energies in relation to system
    %          initialize:
    
    
    
    properties
        
        Points    % mesh points
        Imag_Infinitesimal % imaginary infinitesimal - used for Points_i
        
        Delta     % distance of mesh points (if grid is not uniform distributed, Delta gives an average)
        N         % Number of mesh points
        Limits    % range of mesh
        Uniform   % checks wheater Grid is uniform distributed
    end
    
    
    methods
        
        
        function obj = Grid(grid_points, imag_part)
            if nargin < 2 || isempty(imag_part), imag_part = 10^-5; end
            obj.Points = grid_points(:);
            difference = diff(grid_points);
            
            obj.Limits = [min(grid_points), max(grid_points)];
            obj.N = numel(grid_points);
            obj.Imag_Infinitesimal = imag_part;
            obj.Delta = min(difference);
            if var(difference) < 10^-18
                obj.Uniform = true;
            else
                obj.Uniform = false;
            end
          
        end
          
        
        function obj = f_add(obj, x)
            %TODO: not works for vector input - rethink insertion
            %check if in limits
            if any(x <= obj.Limits(1) )

                obj.Limits(1) = min(x);
            end
            if any(x >= obj.Limits(2))
                obj.Limits(2) = max(x);
                %TODO: change L for vector input
                L = obj.N + 1;
            else
                L = find(obj.Points > x, 1);
            end
            obj.Points = [obj.Points(1:L), x,obj.Points((L+1):end)];
            obj.N = obj.N + numel(x);
            obj.Points_i = obj.Points + obj.Imag_Infinitesimal;
            
            
            % not ready for vector input!
            obj.Delta = (obj.Delta * (obj.N-2) - (obj.Points(L+1)-obj.Points(L-1)) + sum(diff(obj.Points((L-1):(L+1)))))/(obj.N-1);
            %obj.Delta = mean(diff(obj.Points)); % in order to check
            
            
        end
        
        function points_i = Points_i(obj)
            points_i = obj.Points + 1i* obj.Imag_Infinitesimal;
        end
        
        function index = f_find_index(obj,x)
            %finds the index of the grid next to x
            
            if numel(x) > 1
                [gr,xx]= meshgrid(obj.Points,x);
            else
                gr = obj.Points;
                xx = x;
            end
                
            
            [~,index] = min(abs(gr-xx),[],1);
        end
        
        
        function obj = f_refine_by_new_delta(obj, delta)
            if obj.Uniform == true
                interval = diff(obj.Limits);
                num = interval/delta;
                new_points = linspace(obj.Limits(1), obj.Limits(2), num);
                obj = Grid(new_points, obj.Imag_Infinitesimal);
            
            else
                error('grid is not uniform')
            end
        end

        function obj = f_refine(obj, x, refine_number)
            %function to add grid points to refine grid for special values
            % search for position x and refine grid for 10 nodes
            if x <= obj.Limits(1) || x >= obj.Limits(2)
                
            end
            
            warning('not ready')
            
        end
        
    end
    
end
        