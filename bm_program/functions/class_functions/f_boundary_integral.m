function [integ, parameters] = f_boundary_integral(x,y, left, mes_flag)
if nargin < 4 || isempty(mes_flag)
    mes_flag = false;
end
% approximates the integral from -inf to inf outisde the boundaries of the function y with
% the rational function a/(x-b)^2
% the integral is thus +-a/(x_boundary-b)
% left indicates the integral (-inf, left_boundary), not left: (right_boundary, inf)

if left
    x_boundary = x(1);
else
    x_boundary = x(end);
end


slope_abs = mean(diff(abs(y))./diff(x));
if sign(round(slope_abs*10^7)) .* (left-1/2) <= 0
    %TODO: Flag to print warning (mes_flag)
    if mes_flag
        warning('function is not going to zero for |x|-> inf, no integral approximation to inf')
    end
    integ = 0;
    parameters = [0,0];
else
    
    %NOTE: check if y has same argument for most of its values (otherwise difficult)
    
    
    slope = mean(diff(y)./diff(x));
    b = mean( 2*y./slope + x);
    a = mean((x-b).^2.*y);
    
    integ = -(left-1/2)*2*a./(x_boundary-b);
    parameters = [a;b];
    
%     a./(x-b).^2
%     y
end