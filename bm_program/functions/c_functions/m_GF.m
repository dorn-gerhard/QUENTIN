function [GR, GK] = m_GF(points, i0p, E_p, E_n, E_q, C1, C2, C3, C4, rho)


[E_pp, E_nn] = meshgrid(E_p, E_n);
[E_mm,E_qq] = meshgrid(E_n, E_q);
GR = complex(zeros(numel(points),1));
GK = complex(zeros(numel(points),1));
for k = 1:numel(points)
    GR(k) = trace( rho * (1./(points(k) - (E_pp - E_nn)+i0p*1i).*C1) * C2 + rho * C3 * (1./(points(k) - (E_mm - E_qq) + i0p*1i) .* C4));
    GK(k) = trace( -2i*i0p * rho * (1./( (points(k) - (E_pp - E_nn)).^2 + i0p^2) .* C1) * C2  + 2i*i0p * rho * C3 * (1./( (points(k) - (E_mm - E_qq)).^2 + i0p^2) .* C4) );
end

% 
% GR = mmat(mmat(rho, bsxfun(@times, (1./(permute(points,[2,3,1]) - (E_p' - E_n) + i0p*1i)),C1)), C2) + ...
%     mmat(mmat(rho, C3), bsxfun(@times, (1./(permute(points, [2,3,1]) - (E_n' - E_q) + i0p*1i)), C4));
% n = size(GR,3);
% [GR_real, GR_imag, GK_real, GK_imag] = deal(zeros(n,1));
% for k = 1:n
%     GR_real(k,1) = real(trace(GR(:,:,k)));
%     GR_imag(k,1) = imag(trace(GR(:,:,k)));
% end
% 
% GK = -2i*i0p*   mmat(rho, mmat(bsxfun(@times, (1./((permute(points,[2,3,1]) - (E_p' - E_n)).^2 + i0p^2)), C1), C2)) + ...
%     2i*i0p*   mmat(rho, mmat(C3, bsxfun(@times, (1./((permute(points,[2,3,1]) - (E_n' - E_q)).^2 + i0p^2)), C4)));
% 
% 
% for k = 1:n
%     GK_real(k,1) = real(trace(GK(:,:,k)));
%     GK_imag(k,1) = imag(trace(GK(:,:,k)));
% end
% 
% 
% 
% end
% 
% function C = mmat(A,B)
% 
% na = ndims(A);
% nb = ndims(B);
% C = zeros(size(A,1), size(B,2), max(size(A,3), size(B,3)));
% if na == nb
%     %equal size case
%     if na == 2
%         %matrix multiplication
%         C = A*B;
%     elseif na == 3
%         n = size(A,3);
%         for k = 1:n
%             C(:,:,k) = A(:,:,k) * B(:,:,k);
%         end
%     end
% else
%     if na == 3 && nb == 2
%         n = size(A,3);
%         for k = 1:n
%             C(:,:,k) = A(:,:,k) * B(:,:);
%         end
%     elseif na == 2 && nb == 3
%         n = size(B,3);
%         
%         for k = 1:n
%             C(:,:,k) = A(:,:) * B(:,:,k);
%         end
%     end
% end
% 
% end
% 
