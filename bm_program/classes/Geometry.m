classdef Geometry
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N_site
        Bonds        
        description
        N_bonds
        Bond_index
        Coordinates = []
    end
    
    methods
        %constructor
        function obj = Geometry(N_site, type, bonds, coordinates)
            if nargin < 1 || isempty(N_site), error('number of sites has to be specified'); end
            if nargin < 2 || isempty(type) || ~ischar(type), error('type has to be specified: p..periodic, c..cyclic, f..full'); end
            if nargin < 4 || isempty(coordinates), coordinates = []; end
            
            switch lower(type(1))
                case 'f' %full
                    if nargin < 3, error('full bonding structure must be given'); end
                    if ~ismatrix(bonds) || size(bonds,1) ~= size(bonds,2) || size(bonds,1) ~= N_site, error('bonds must be square matrix (N_site x N_site)'); end
                    obj.Bonds = bonds;
                    obj.description = 'full bonding structure';
                    obj.Bond_index = reshape(1:N_site^2,N_site, N_site);
                case 'p'  %periodic, cyclic%
                    obj.Bonds = diag(ones(N_site-1,1),1) + diag(ones(N_site-1,1),-1);
                    obj.Bonds(1,end) = 1;
                    obj.Bonds(end,1) = 1;
                    obj.description = 'periodic cycle';
                    obj.Bond_index = diag(1:(N_site-1),1) + diag(-1:-1:-(N_site-1),-1);
                    obj.Bond_index(1,end) = N_site;
                    obj.Bond_index(end,1) = -N_site;
                case 'l' % chain structure
                    obj.Bonds = diag(ones(N_site-1,1),1) + diag(ones(N_site-1,1),-1);
                    obj.description = 'linear chain';
                    obj.Bond_index = diag(1:(N_site-1),1) + diag(-1:-1:-(N_site-1),-1);
                otherwise
                    error('not recognized geometry input, available: f (full), p (periodic), l (chain)')
                                            
            end
            obj.N_bonds = sum(obj.Bonds(:) ~= 0);
            obj.N_site = N_site;
            obj.Coordinates = coordinates;
        end
        
        function next_nearest_neighbour = nnn(obj,distance)
            %calculate distances
            if nargin < 2 || isempty(distance) || ~isnumeric(distance) || numel(distance) ~= 1, distance = 2; end
            next_nearest_neighbour = obj.Bonds^distance;
        end
        
        function s_geometry(obj, bond_label)
            if nargin < 2 || isempty(bond_label), bond_label = obj.Bonds; end
            
            if isempty(obj.Coordinates), error('no coordinates defined yet.'); end
            %try to find a way to characterize geoemtry, first of all check if it's a planar graph
            gplot(obj.Bonds, obj.Coordinates,'*--')
            axis equal
            for k = 1:obj.N_site
                text(obj.Coordinates(k,1),obj.Coordinates(k,2),num2str(k), 'verticalalignment', 'bottom', 'horizontalalignment', 'center')
                L = true(1,obj.N_site);
                L(1:k) = false;
                for j = find((bond_label(k,:)~= 0) & L)
                    A = obj.Coordinates(k,:);
                    B = obj.Coordinates(j,:);
                    C = (A+B)/2;
                    
                    text(C(1), C(2), num2str(bond_label(k,j),2), 'color', 'r')
                end
            end
            title('hopping direction t_{ij} for i<j')
            
            plot_settings()
            
        end
        
        function f_planar_graph_coordinates_from_bonds()
            disp('not realized yet - not trivial problem');
        end
        
    end
    
end

