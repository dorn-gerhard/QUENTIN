

A = [rand(1,3), zeros(1,2), zeros(1,2); zeros(2,3), rand(2), zeros(2,2); zeros(2,3), zeros(2,2), rand(2)];
A = A + 1i*rand(size(A))
Er = [1,1,1,2,2];

AA = kron(conj(A),A)

n = 50;
bll = rand(10,1);
bll2 = ceil(bll/sum(bll)*n);
sum(bll2) - n

block_size = [3;2;2]
additive = [0; cumsum(block_size(1:end-1).^2)];

Ind = arrayfun(@(x,y) reshape(y+(1:x.^2), x,x),block_size, additive, 'Uniformoutput', false);
index = blkdiag(Ind{:})



L = index >= 1;

AF = AA(:,L(:))

%kann man dieses Resultat auch anders erreichen?
% JA => blockwise


for k = 1:numel(block_size)
    blk{k} = (1:block_size(k)) + sum(block_size(1:k))-block_size(k); 
    AI{k} = kron(conj(A(:,blk{k})), A(:,blk{k}));
end
result = cell2mat(AI(1:end));


B = AF - cell2mat(AI(1:end));
sum(B(:))


%Minimierungsaufgabe von Blockgrößen


Tab2 = Table(En, 'NSE')
En = condense_energy_blocks(En, 0.4)

Tab2 = Table(En, 'NSA')
mean(Tab2.Block_Size)
bar(Tab2.Block_Size)

Qtab = get_Qtab(En)

Qtab.s_print

bar(Tab2.Block_Size('A', {-300,-95}))
Tab2.Average_Energies

