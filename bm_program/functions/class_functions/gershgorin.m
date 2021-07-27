function gershgorin(A, plot_ew)
if nargin < 2 || isempty(plot_ew), plot_ew = false; end



t = linspace(0,2*pi,400);



for k = 1:length(A)
    r = min(sum(  abs(A(k,[1:k-1, k+1:end]))), sum(  abs(A([1:k-1, k+1:end, k]))));
    plot(r*cos(t) + real(A(k,k)), r*sin(t) + imag(A(k,k)),'--')
    axis equal
    hold on
end


if plot_ew
    ew = eig(A);
    plot(real(ew), imag(ew), 'xk')
end
    