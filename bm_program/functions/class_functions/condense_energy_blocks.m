function En = condense_energy_blocks(En, Tolerance)
% Gerhard Dorn, 17.11.2018
% 
% finds out if blocks should be fused according to tolerance
% NOTE: adopted for spinless systems and systems without perserved spin
% particle / spin block structure (preserved quantity) is stored in Energy.Structure!!!

if En.Basis.N_spin > 1 %distinguish between different block structure depending
    %on conservation of particle number (N) and spin (S)
    if ismember('N', En.Structure) && ismember ('S', En.Structure) 
        Qtab = Table(En,'NS');
    elseif ismember('N', En.Structure)
        Qtab = Table(En, 'N');
    elseif ismember('S', En.Structure)
        Qtab = Table(En, 'S');
    else
        Qtab = Table(En);
        
    end
    
    
else %spinless
    if ismember('N', En.Structure)
        Qtab = Table(En, 'N');
    else
        Qtab = Table(En);
    end
end

plot_flag = false;


New_block = zeros(En.Numel_En,1);
Average_Energy = zeros(En.Numel_En,1);
Gap_Energy = zeros(En.Numel_En,1);
for index = 1:Qtab.Numel_Blocks
    if ismember('N', En.Structure) && ismember ('S', En.Structure)
        N_part = Qtab.N_part(index);
        Spin = Qtab.Spin(index);
        Ener = En.Energies('NS', N_part, Spin);
        en_index = En.Index('NS', N_part, Spin);
    elseif ismember ('S', En.Structure)
        Spin = Qtab.Spin(index);
        Ener = En.Energies('S', Spin);
        en_index = En.Index('S', Spin);
    elseif ismember('N', En.Structure)
        N_part = Qtab.N_part(index);
        Ener = En.Energies('N', N_part);
        en_index = En.Index('N', N_part);
    else
        Ener = En.Energies;
        en_index = En.Index;
    end
    
    LL = abs([0;diff(Ener(:))]) < Tolerance;
    
    if index == 1
        New_block_index = cumsum(~LL)+1;
    else
        New_block_index = cumsum(~LL)+1 + New_block(min(en_index)-1);
    end
    
    New_block(en_index) = New_block_index;
    
    block_size =  histcounts(New_block_index, (min(New_block_index):1+max(New_block_index))-0.5).';
    En_mean = cellfun(@(x) mean(x), mat2cell(Ener, block_size, 1));
    Gap = cellfun(@(x) max(x) - min(x), mat2cell(Ener, block_size,1));
    
    Average_Energy(en_index) = cell2mat(arrayfun(@(x,y) ones(x,1) * y,...
        block_size, En_mean, 'uniformoutput', false));
    Gap_Energy(en_index) = cell2mat(arrayfun(@(x,y) ones(x,1) * y,...
        block_size, Gap, 'uniformoutput', false));
    if plot_flag
        errorbar((N_part + 0.08*Spin)*ones(size(en_index)), Average_Energy(en_index),Gap_Energy(en_index)/2, 'ob')
        hold on 
        plot((N_part + 0.08*Spin)*ones(size(en_index)), Ener, 'xr')
        grid on
    end
    
end

if En.Categ.List.isKey('Average_Energies') || En.Categ.List.isKey('A') 
    En.Categ = En.Categ.f_remove({'Average_Energies'});
   
end
if En.Categ.List.isKey('Gap_Energies') || En.Categ.List.isKey('G')
     En.Categ = En.Categ.f_remove({'Gap_Energies'});
end
En.Categ = En.Categ.f_append({'Average_Energies', 'A'}, Average_Energy);
En.Categ = En.Categ.f_append({'Gap_Energies', 'G'}, Gap_Energy);

En.Structure = En.Categ.f_shortnames;


end



