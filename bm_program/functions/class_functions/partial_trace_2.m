function sigma = partial_trace_2(rho, dimension, sub1_or_2)
% e.g. rho  .. 6x6
% dimension ..  (2,3)
% sub1_or_2 .. 1
% result: 3x3 matrix

    superblock = mat2cell(rho, repmat(dimension(2),dimension(1),1), repmat(dimension(2),dimension(1),1));
if sub1_or_2 == 2    
    sigma = cellfun(@(x) trace(x), superblock);
elseif sub1_or_2 == 1
    diag_blocks = arrayfun(@(x) superblock(x,x), 1:dimension(1), 'uniformoutput', true);
    sigma = sum(reshape(cell2mat(diag_blocks), dimension(2), dimension(2), dimension(1)),3);
end


