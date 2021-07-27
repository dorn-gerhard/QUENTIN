function [Qtab] = get_Qtab(En)
if En.Basis.N_spin == 2,
    spin_used = 'S';
else
    spin_used = [];
end

Qtab = Table(En,['N', spin_used]);

Ground_en = zeros(Qtab.Numel_Blocks,1);
for index = 1:numel(Qtab.Index);
    block_index = Qtab.Index(index);
    N_el = Qtab.N_part(block_index);
    Spin = Qtab.Spin(block_index);
    Ground_en(index) = min(En.Energies('NS', N_el, Spin));
end
Qtab = Qtab.f_add_categ({'Groundenergy', 'G'}, Ground_en);