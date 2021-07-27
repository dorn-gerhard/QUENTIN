function g = cauchy(om,f, delta)
if nargin < 3 || isempty(delta), delta = 10^-5; end

%  y is a matrix
% g(x) = P \int_-inf^inf f(omega)/ (omega - x) d\omega
% = 2 \delta f'(x) + \int_-inf^{x-delta} f(\omega) / (omega - x) dx + 
%   \int_{x+delta}^inf f(\omega) / ( \omega - x) d\omega 


% split x in tranches if meshgrid becomes to large
om = om(:);
f = f(:);
x = om.';

% for k = 1:m
% max_size = 1000 (e.g.)
% or define Index set (last set has not max_size elements!
% x = w((1:max_size) + max_size * (k-1))



%create integration weights!

N = numel(om);

if mod(N,2) == 1
    x_a = om(1:2:end-1);
    x_b = om(2:2:end);
    x_c = om(3:2:end);
else
    x_a = om(1:2:end-2);
    x_b = om(2:2:end-1);
    x_c = om(3:2:end-1);
end


w_a = -(x_a - x_c) .*(2*x_a - 3*x_b + x_c)./(6*(x_a-x_b));
w_b = (x_a - x_c).^3./(6*(x_b - x_a) .*(x_b - x_c));
w_c = (x_a - x_c) .*(x_a - 3*x_b + 2*x_c)./(6*(x_b - x_c));
N_interval = numel(w_a);  %number of intervals

w = zeros(size(om));

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

% create integration array

[xx,oo] = meshgrid (x,om);
ind_x = 1:numel(x);

L = xx == oo; % exclude those integrals
% identify which w- value corresponds to the pole?
[pole_index,~] = ind2sub(size(L),find(L));
interval_number = floor(pole_index/2).'; %interval number
is_interval_boundary = logical(mod(pole_index,2)'); % 0 in between, 1 on boundary

if any(interval_number > N_interval)
    L_bad_interval = interval_number > N_interval;
    interval_number(L_bad_interval) =  [];
    is_interval_boundary(L_bad_interval) = [];
    %pole_index(L_bad_interval) = [];
    treat_last_interval_flag = true;
else
    treat_last_interval_flag = false;
end

ww = repmat(w,1,numel(x));


% remove integral weights for all first intervals (boundary and not boundary)
int_num_1 = interval_number(interval_number >= 1);
ind_x_1 = ind_x(interval_number >= 1);
subtract = [w_a(int_num_1).';w_b(int_num_1).'; w_c(int_num_1).'];
index = sub2ind(size(ww), (int_num_1-1)*2 + (1:3)', repmat(ind_x_1,3,1));
ww(index) = ww(index) - subtract;

% remove integral weights for next interval
int_num_2 = interval_number(interval_number + 1 <= N_interval & is_interval_boundary); %check that number of intervals is not exceeded
ind_x_2 = ind_x(is_interval_boundary & interval_number + 1 <= N_interval);
subtract_2 = [w_a(int_num_2+1).';w_b(int_num_2+1).'; w_c(int_num_2+1).'];
index_2 = sub2ind(size(ww), int_num_2*2 + (1:3)', repmat(ind_x_2,3,1));
ww(index_2) = ww(index_2) - subtract_2;

% add lower trapezoidal rule for boundary elements
int_num_3 = interval_number(interval_number >= 1 & is_interval_boundary);
ind_x_3 = ind_x(is_interval_boundary & interval_number >= 1);
row_add_lower = (int_num_3-1)*2+(1:2)';
trapz_weight = repmat(diff(om(row_add_lower))/2,2,1);
index_3 = sub2ind(size(ww), row_add_lower, repmat(ind_x_3,2,1));

ww(index_3) = ww(index_3) + trapz_weight;

% add higher trapezoidal rule for boundary elements
int_num_4 = int_num_2;
ind_x_4 = ind_x_2;
row_add_higher = int_num_4 * 2 + (2:3)';
trapz_weight = repmat(diff(om(row_add_higher)/2),2,1);
index_4 = sub2ind(size(ww), row_add_higher, repmat(ind_x_4,2,1));

ww(index_4)= ww(index_4) + trapz_weight;


%interpolate f(\om) on points: pole +- delta
%create interpolation polynomial using points x-1, x, x+1 using Lagrange polynomials:
L_pole_not_edge = pole_index - 1 >= 1 & pole_index +1 <= N;
    
a = om(pole_index(L_pole_not_edge)-1).';
b = om(pole_index(L_pole_not_edge)).';
c = om(pole_index(L_pole_not_edge)+1).';
fa = f(pole_index(L_pole_not_edge)-1).';
fb = f(pole_index(L_pole_not_edge)).';
fc = f(pole_index(L_pole_not_edge)+1).';

f_plus_delta = delta * (b+ delta - c)./((a-b).*(a-c)) .* fa + (b+delta - a).*(b+delta - c)./((b-a).*(b-c)) .* fb + (b+delta - a)*delta ./((c-a).* (c-b)) .*fc;
f_minus_delta = -delta * (b- delta - c)./((a-b).*(a-c)) .* fa + (b-delta - a).*(b-delta - c)./((b-a).*(b-c)) .* fb + (b-delta - a)*-delta ./((c-a).* (c-b)) .*fc;


% add trapzrule for those 4 points (f_{x-1}, f_{x-\delta}, f_{x+\delta}, f_{x+1}) and add trapezoidal rule
ww_plus = zeros(4,numel(x));
ww_plus(:,L_pole_not_edge) = [repmat(b-delta-a,2,1); repmat(c-b-delta,2,1)]/2;
%integrand:

ff = f./(om-x);%;.*sign(om-x);
ff(isinf(ff)) = 0;
ff(isnan(ff)) = 0;

%if integrand ff is not quite zero at the edge of the grid, extrapolate with 1/x
% add integration of limits to -inf +inf by
% Strategy: take last grid points and fit a 1/x function which can be integrated
% analytically.
L_left_boundary = 1:10;
L_right_boundary = numel(om) + (-9:0);
left_integrand =  ff(L_left_boundary,:);
right_integrand = ff(L_right_boundary,:);

mes_flag = true;
left = true;

a_left = mean(left_integrand.*(om(L_left_boundary) - x).^2,1);
left_boundary = om(1);
Int_left = -a_left./(left_boundary - x);

a_right = mean(right_integrand.*(om(L_right_boundary) - x).^2,1);
right_boundary = om(end);
Int_right = a_right./ (right_boundary - x);

%[left_integ,  left_parameters] = f_boundary_integral(om(L_left_boundary), left_integrand, left, mes_flag);
%[right_integ, right_parameters] = f_boundary_integral(om(L_right_boundary), right_integrand, ~left, mes_flag);
% 
% plot(om, f)
% grid on 
% hold on
% 
% t = linspace(-150,-110,1000);
% plot(om, ff(:,100))
% hold on
% 
% plot(t, a_left(100)* 1./(t- x(100)));
% 
% z = linspace(110, 150, 1000);
% plot(z, a_right(100)*1./(z-x(100)));
% 
% left_boundary - x


ff_plus = zeros(4, numel(x));
ff_plus(:,L_pole_not_edge) = [fa./(a-b); f_minus_delta./(-delta); f_plus_delta./(delta); fc./(c-b)];

ww_gesamt = [ww;ww_plus];
ff_total = [ff; ff_plus];





% calculate derivative of f at poles and add to integral
% use derivative of Lagrange polynom from a,b,c:
f_prime = zeros(1,numel(x));
f_prime(L_pole_not_edge) = (b-c)./((a-b).*(a-c)) .* fa + (1./(b-a) + 1./(b-c)).*fb + (b-a)./((c-a).*(c-b)).*fc;






g =sum(ff_total.*ww_gesamt)+2 * delta * f_prime + Int_left + Int_right;


