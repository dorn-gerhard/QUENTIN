classdef Test
    
    properties
        
        
       A 
       B 
       C
    end
    
    methods
        function obj = Test(struc)
            struct2vars(struc)
            if ~isfield(struc,'A'), A = 0; end
            if ~isfield(struc,'B'), B = 3; end
            if ~isfield(struc,'C'), C = 2; end
            obj.A = A;
            obj.B = B;
            obj.C = C;
        
            
            
            
        end
    end
end