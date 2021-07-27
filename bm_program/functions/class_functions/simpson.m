function I = simpson(x,y)


N = numel(x);

if mod(N,2) == 1
    x_a = x(1:2:end-1);
    x_b = x(2:2:end);
    x_c = x(3:2:end);
else
    x_a = x(1:2:end-2);
    x_b = x(2:2:end-1);
    x_c = x(3:2:end-1);
end


w_a = -(x_a - x_c) .*(2*x_a - 3*x_b + x_c)./(6*(x_a-x_b));
w_b = (x_a - x_c).^3./(6*(x_b - x_a) .*(x_b - x_c));
w_c = (x_a - x_c) .*(x_a - 3*x_b + 2*x_c)./(6*(x_b - x_c));


w = zeros(size(x));

if mod(N,2) == 1
    w(1:2:end-1) = w_a;
    w(2:2:end) = w_b;
    w(3:2:end) = w(3:2:end) + w_c;
else
    w(1:2:end-2) = w_a;
    w(2:2:end-1) = w_b;
    w(3:2:end-1) = w(3:2:end-1) + w_c;
    w(end-1) = w(end-1) + 1/2*(x(end) - x(end-1));
    w(end) = 1/2*(x(end) - x(end-1));
end

I =sum(y.*w);


