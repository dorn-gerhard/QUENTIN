function g = cauchy_trapz(om,f, delta)
if nargin < 3 || isempty(delta), delta = 10^-5; end
if delta > min(diff(om))
    warning('chosen delta is smaller than smallest grid interval - use smallest grid interval*3/4 as delta!')
    delta = min(diff(om))*3/4;
end
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

x_a = om(1:end-1);
x_b = om(2:end);

w_a = (x_b - x_a)/2;
w_b = (x_b - x_a)/2;
N_interval = numel(w_a);  %number of intervals
w = zeros(size(om));

ww_a = [repmat(w_a,1,numel(x)); zeros(1,numel(x))];
ww_b = [zeros(1,numel(x)); repmat(w_b,1,numel(x)) ];

w(1:end-1) = w_a;
w(2:end) = w(2:end) + w_b;
 
ww = repmat(w,1,numel(x));

% %%
% %delete integration weights for intervals within the radius:
% radius = 1;
% 
% %delete all integral weights of the intervals in between for a symmetric integration of
% %(pole - radius) to (pole + radius) leaving out the delta area 
% %delete_all_those_data_points +-1 
% L = abs(x-om) < 1;
% %find the number of intervals, that have a distance more than the radius
% 
% ww_a(L) = 0;
% ww_b(L) = 0;
% lower_int_index = sum(x-om > radius,1);
% upper_int_index = sum(x-om +radius > 0,1)+1;
% LL_lower = (lower_int_index > 0) & (lower_int_index <= N);
% LL_upper = (upper_int_index > 0) & (upper_int_index <= N);
% LL = LL_lower & LL_upper;
% ww_a(sub2ind([N,N], lower_int_index(LL_lower), find(LL_lower))) = 0;
% ww_b(sub2ind([N,N], upper_int_index(LL_upper), find(LL_upper))) = 0;
% 
% 
% ww = ww_a + ww_b;
% 
% [xx,oo] = meshgrid (x,om);
% 
% ff = f./(om-x);%;.*sign(om-x);
% ff(L) = 0;
% 
% % 2) integrate from interval boundary to pole +- 1 via trapezoidal rule
% 
% % lower point: 
% 
% L_pole = xx == oo; % exclude those integrals
% % identify which w- value corresponds to the pole?
% [pole_index,~] = ind2sub(size(L_pole),find(L_pole));
% L_pole_not_edge = pole_index > 1 & pole_index < N;
% 
% a = om(lower_int_index(LL)).';
% b_low = om(lower_int_index(LL)+1).';
% b_inter_low = om(pole_index(LL)).' - radius;
% fa = f(lower_int_index(LL)).';
% fb_low = f(lower_int_index(LL)+1).';
% f_lower = (b_inter_low  - b_low) ./(a-b_low) .* fa + (b_inter_low - a)./(b_low-a).* fb_low ;
% 
% 
% c = om(upper_int_index(LL)).';
% b_inter_upp = om(pole_index(LL)).' + radius;
% b_upp = om(upper_int_index(LL)-1).';
% fb_upp = f(upper_int_index(LL)-1).';
% fc = f(upper_int_index(LL)).';
% 
% f_upper = (b_inter_upp - b_upp) ./ (c - b_upp) .* fc + (b_inter_upp - c)./(b_upp - c) .* fb_upp;
% 
% ff_add = zeros(4,N);
% ff_add(1:2,LL) = [fa./(a-b_inter_low-radius); f_lower./(-radius)];
% ff_add(3:4,LL) = [f_upper./radius; fc./(c-b_inter_upp+radius)];
% 
% ww_add = zeros(4,N);
% ww_add(1:2, LL) = repmat(b_low - a, 2,1)/2;
% ww_add(3:4, LL) = repmat(c - b_upp,2,1)/2;
% 
% % 3) integrate symmetrically from delta to radius (f(x) - f(-x))/(x-omega)
% 
% %integrate from b_low=b_upp + delta to b_inter_upp
% n =10;
% for k = 1:numel(LL)
%     low_ind = lower_int_index(LL);
%     upp_ind = upper_int_index(LL);
%     ivec = [b_inter_low(k); om(low_ind(k):upp_ind(k)-1); b_inter_upp(k)]
%     pole = b_upp(k)-radius
% 
% 
% vec = delta + (radius-delta)/n.*(0:n)'; %add b_low
% 
% %interp1(vec, 
% 
% %%
% 4) integrate pole via derivative

% 5) integrate to +- infinity

% create integration array

[xx,oo] = meshgrid (x,om);

ff = f./(om-x);%;.*sign(om-x);
ff(isnan(ff)) = 0;
ff(isinf(ff)) = 0;

% 2) integrate from interval boundary to pole +- 1 via trapezoidal rule

% lower point: 

L_pole = xx == oo; % exclude those integrals
% identify which w- value corresponds to the pole?
[pole_index,~] = ind2sub(size(L_pole),find(L_pole));
L_pole_not_edge = pole_index > 1 & pole_index < N;


ind_x = 1:numel(x);


interval_number = floor(pole_index).'; %interval number


if any(interval_number > N_interval)
    L_bad_interval = interval_number > N_interval;
    interval_number(L_bad_interval) =  [];
    %pole_index(L_bad_interval) = [];
    treat_last_interval_flag = true;
else
    treat_last_interval_flag = false;
end



%remove integration weights around pole:
ww = ww - diag(w_a,1) - diag([0;w_b])  - diag([w_a;0],0) - diag(w_b,-1);

L_pole_not_edge = pole_index > 1 & pole_index < N;

a = om(pole_index(L_pole_not_edge)-1).';
b = om(pole_index(L_pole_not_edge)).';
c = om(pole_index(L_pole_not_edge)+1).';
fa = f(pole_index(L_pole_not_edge)-1).';
fb = f(pole_index(L_pole_not_edge)).';
fc = f(pole_index(L_pole_not_edge)+1).';


n =20;
f_plus_delta = (b+delta - b) .* (b+ delta - c)./((a-b).*(a-c)) .* fa + (b+delta - a).*(b+delta - c)./((b-a).*(b-c)) .* fb + (b+delta - a).*(b+delta-b) ./((c-a).* (c-b)) .*fc;
f_minus_delta = (b-delta  - b) .* (b- delta - c)./((a-b).*(a-c)) .* fa + (b-delta - a).*(b-delta - c)./((b-a).*(b-c)) .* fb + (b-delta - a).*(b-delta - b)./((c-a).* (c-b)) .*fc;

a_vec = a + (b-delta-a)/n.*(0:n)';
c_vec = b+delta + (c-(b+delta))/n.*(0:n)';
f_minus_delta_vec = (a_vec  - b) .* (b- delta - c)./((a-b).*(a-c)) .* fa + (a_vec - a).*(b-delta - c)./((b-a).*(b-c)) .* fb + (a_vec - a).*(b-delta - b)./((c-a).* (c-b)) .*fc;
f_plus_delta_vec = (c_vec - b) .* (c_vec - c)./((a-b).*(a-c)) .* fa + (c_vec - a).*(c_vec - c)./((b-a).*(b-c)) .* fb + (c_vec - a).*(c_vec - b) ./((c-a).* (c-b)) .*fc;

w_a_delta = repmat((b-delta-a)/n,n+1,1);
w_a_delta([1,end],:) = w_a_delta([1,end],:)/2;

w_c_delta = repmat((c-(b+delta))/n, n+1, 1);
w_c_delta([1,end],:) = w_c_delta([1,end],:)/2;

ff_minus = f_minus_delta_vec./(a_vec - b);
ff_plus = f_plus_delta_vec./ (c_vec -b );

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


%ff_plus = zeros(4, numel(x));
%ff_plus(:,L_pole_not_edge) = [fa./(a-b); f_minus_delta./(-delta); f_plus_delta./delta; fc./(c-b)];

ff_add = zeros(2*n+2, numel(x));
ff_add(:,L_pole_not_edge) = [ff_minus;ff_plus];
ww_plus = zeros(2*n+2, numel(x));
ww_plus(:,L_pole_not_edge) = [w_a_delta; w_c_delta];

ww_gesamt = [ww;ww_plus];
ff_total = [ff; ff_add];





% calculate derivative of f at poles and add to integral
% use derivative of Lagrange polynom from a,b,c:
f_prime = zeros(1,numel(x));
f_prime(L_pole_not_edge) = (b-c)./((a-b).*(a-c)) .* fa + (1./(b-a) + 1./(b-c)).*fb + (b-a)./((c-a).*(c-b)).*fc;






g = sum(ff_total.*ww_gesamt)+2 * delta * f_prime; %+ Int_left + Int_right;


