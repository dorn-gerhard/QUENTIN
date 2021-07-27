classdef Sigma_grandcanonical < Sigma
    
    %author: Gerhard Dorn
    %date: october 20
    
    
    properties
        Beta
        Mu
    end
    
    
    methods
        function obj = Sigma_grandcanonical(Mu_inp, Beta_inp, Ener_inp, Tab_inp, Energy_Cut_inp, Energy_Cut_2_inp, sigma_tolerance)
            
            %% Pre Initialization %%
            if nargin < 7 || isempty(sigma_tolerance), sigma_tolerance = 10^-7; end
            if nargin < 5 || isempty(Energy_Cut_inp), Energy_Cut_inp = 10^10; end
            if nargin < 6 || isempty(Energy_Cut_2_inp), Energy_Cut_2_inp = Energy_Cut_inp; end
            if nargin < 4 || isempty(Tab_inp), Tab_inp = Table(Ener_inp, 'NE'); end
            if Tab_inp.Numel_Datasets ~= Ener_inp.Numel_En
                warning('Table does not describe full range of Energies - therefore new full table is generated')
            end
            %new strategy!!!
            %create diagonal and turn this into the block structure with is
            %defined via Tab
            %Then cut those blocks that are below the tolerance
            
            rho_diagonal = exp(-Beta_inp .*(Ener_inp.Energies - min(Ener_inp.Energies) - Mu_inp*Ener_inp.N_part));
            rho_diagonal(abs(rho_diagonal) < sigma_tolerance) = 0;
            rho_diagonal = sparse(diag(rho_diagonal)./sum(rho_diagonal));
            
            blocks = Tab_inp.f_block(rho_diagonal);
            
            % delete blocks
            block_index = 1:numel(blocks);
            blocks_to_delete = cellfun(@(x) sum(diag(full(x)))< sigma_tolerance/100, blocks);
            block_index(blocks_to_delete) = [];
            blocks(blocks_to_delete)=[];
            new_tab = Tab_inp.f_cut_filter('I', block_index);
            
            
            
            
            %% Object Initialization %%
            obj = obj@Sigma(Ener_inp, blocks, new_tab, Energy_Cut_inp,Energy_Cut_2_inp);
            
            
            %% Post Initialization %%
            obj.Beta = Beta_inp;
            obj.Mu = Mu_inp;
        
        end
        
    end
    
    
end