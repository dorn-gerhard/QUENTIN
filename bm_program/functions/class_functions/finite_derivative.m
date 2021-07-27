function dx = finite_derivative(f, dimen, h)
if nargin < 2 || isempty(dimen), dimen = 1; end

% use 5 point stencil
%0,1,2,3,4: 
%f'_x = (-25*f[i+0]+48*f[i+1]-36*f[i+2]+16*f[i+3]-3*f[i+4])/(12*h)
%-1, 0, 1, 2, 3
%f'_x = (-3*f[i-1]-10*f[i+0]+18*f[i+1]-6*f[i+2]+1*f[i+3])/(12*h)
%-2 -1 0 1 2
% f_x = (1*f[i-2]-8*f[i-1]+0*f[i+0]+8*f[i+1]-1*f[i+2])/(12*h)
%-3 -2 -1 0 1
%f_x = (-1*f[i-3]+6*f[i-2]-18*f[i-1]+10*f[i+0]+3*f[i+1])/(12*h)
%-4 -3 -2 -1 0
%f_x = (3*f[i-4]-16*f[i-3]+36*f[i-2]-48*f[i-1]+25*f[i+0])/(12*h)

dims = size(f);



if dimen == 2
    f = f.';
    dx = zeros(size(f));
end

if numel(h) == 1 % equidistant grid points
    %dx(1,:) = [-25, 48, -36, 16, -3] * f(1:5,:) / ( 12*h);
    if size(f,1) < 5
        dx(1:end-1,:) = diff(f,1,1);
    else
        dx(1,:) = [-3, 4, 1] * f(1:3,:) / (2*h);
        dx(2,:) = [-3 -10 18 -6 1] * f(1:5,:) / (12*h);
        dx(3:end-2,:) = (1 * f(1:end-4, :) -8 * f(2:end-3,:) + 0 * f(3:end-2,:) + 8 * f(4:end-1,:) -1 * f(5:end,:))/ ( 12 * h);
        dx(end-1,:) = [-1 6 -18 10 3] * f(end-4:end,:) / (12*h);
        %dx(end,:) = [3 -16 36 -48 25] * f(end-4:end,:) / (12*h);
        dx(end,:) = [1 -4 3] * f(end-2:end,:) / (2*h);
    end
else %non equidistant grid points (rescale / reformulate everything)
    error('to be implemented!')
end

if dimen== 2
    dx = dx.';
end

% test = spline((0:numel(f)-1)*h, f);
% 
% figure
% t = linspace(0, numel(f) * h, 200);
% plot(t, ppval(test, t))
% hold on
% plot(t, [diff(ppval(test,t))/(t(2)-t(1)), 0])
% plot((0:numel(f)-1)*h, dx, 'x')






