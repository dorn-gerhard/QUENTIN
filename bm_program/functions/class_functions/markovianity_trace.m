function [nu_m, nu_t, nu_m_t] = markovianity_trace(K, Tab_cut)
%nu_m returns a measure of violation of Markovianity (if it is negative then trace is not preserved
%nu_t returns a measure of violation of trace preservation 
%nu_m_t returns a measure of violation of complete positivity (always positive)
%NOTE: nu_m = nu_m_t + nu_t > 0! nu_m_t > - nu_t
epsilon = 10^-5;

n = size(K,1);

KK = Tab_cut.f_choi(eye(n)+K*epsilon);
[i1,i2,w] = find(KK);
[~,~,rank1] = unique(i1);
[~,~,rank2] = unique(i2);
K2 = sparse(rank1,rank2,w);

%K2 = KK(any(KK,2), any(KK,1));


p = symamd(K2);
K3 = K2(p,p);

%recognizing blocks
%NOTE: add identity matrix for block identification reasons (assume hermiticity), substract in "blocks = ..."
[i2,j2,s2] = find(K3+speye(size(K3)));
%spy(K3)


%NOTE: presort with symamd - then according to the sparse ordering a block on a matrix with a non-zero
%diagonal is found if two diagonal elements occur subsequently, thus the border of a block is found,
%cumsum generates the block number
diagonal_match = i2 == j2;
block_id = cumsum([1;diagonal_match(1:end-1)] == diagonal_match & diagonal_match == 1);
numel_block_elements = hist(block_id, 1:max(block_id)).';
blocks = arrayfun(@(x,y) sparse(i2(x+(1:y))-i2(x+1)+1, j2(x+(1:y)) - j2(x+1) + 1, s2(x+(1:y))) - speye(numel(unique(i2(x+(1:y))))), cumsum([0;numel_block_elements(1:end-1)]), numel_block_elements, 'uniformoutput', false);
%NOTE: control blocks
%AA = blkdiag(blocks{:});

choi_ev = cell2mat(cellfun(@(x) eig(full(x)), blocks, 'uniformoutput', false));
            
choi_sum_ev = sum(choi_ev);
            
choi_pos_ev = sum(abs(choi_ev) - choi_ev);

%can be negative
nu_m_t = (sum(abs(choi_ev)) - Tab_cut.Numel_Datasets)/epsilon;

%violation of positivity normalized by time - always positive
nu_m = choi_pos_ev/epsilon; % equal to taking the trace as reference: (sum(abs(choi_ev)) - choi_sum_ev)/epsilon

%violation of trace:
nu_t = (Tab_cut.Numel_Datasets - choi_sum_ev)/epsilon;



% OLD: 
%ew = eig(full(KK));
%nu_m = (sum(abs(ew)) - Tab_cut.Numel_Datasets)/epsilon;

% NOTE: alternative - very inaccurate!!!

%K2 = KK*KK';
%K2 = K2(any(K2,2), any(K2,1));
%abs(trace(K2^(1/2)) - Tab_cut.Numel_Datasets)/epsilon


%abs(trace(sqrtm(full(K2))) - Tab_cut.Numel_Datasets)/epsilon


