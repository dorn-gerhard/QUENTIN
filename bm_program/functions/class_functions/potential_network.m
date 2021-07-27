
function potential = potential_network(system, contact, pot_bath)

hopping_resistivity = 2; % check which model fits better!


N_bath = size(contact.coupling,1);
if numel(pot_bath) ~= N_bath
    error('Error due to wrong number of elements in pot_bath!')
end



    
N_site = length(system.b);

% create large adjacency matrix

Adj = system.b;

for k = 1:N_bath
    % averaging different bath orbitals
    Adj = [Adj, [sqrt(sum(abs(contact.coupling{k}).^2,2)); zeros(k-1,1)]; ...
           sqrt(sum(abs(contact.coupling{k}).^2,2))', zeros(1,k)];
end



pot_ext = [zeros(size(system.b,1),1);pot_bath(:)];

% Knoten
% unbekannte sind Potentiale an Knoten
Adjac = abs(Adj)> 0 & eye(size(Adj)) == 0;
R = abs(Adj).^hopping_resistivity;
R(~Adjac) = 0;



A = -diag(sum(R,2)) + R; % Knotenregel

A(N_site+(1:N_bath),:) = [];
A(N_site+(1:N_bath), N_site+(1:N_bath)) = eye(N_bath);


potential = A\pot_ext;
end

%Stromrichtungen
%[a,b] = meshgrid(1:length(Adjac));

%R = R.*sign(a-b);

% Maschenregel

% A = [A; zeros(N_bath, size(A, 2)-N_bath), eye(N_bath)]
% 
% A = [A; [A(5,1), 0, 0, - A(5,4), A(5,4), -A(6,1)]]
% 
% U = A\ [zeros(size(A,2),1); Vb/2; -Vb/2; 0]
% %A = [A; zeros(N_bath-1, size(A, 2)-N_bath), 1, -1]
% 
% 
% n = basis.N_site + N_bath;
% 
% off_set = [0, cumsum(n-1:-1:1)]
% f = @(x,y) (off_set(x) + y-x) .* (x < y)
% [yy,xx] = meshgrid(1:n);
% Index = f(xx,yy);
% Index = Index + Index';
% 
% B = zeros(0, nchoosek(n,2))
% 
% for k = 1:n-N_bath
%     neighbours = find(Adjac(k,:) == 1);
%     B(k, Index(k,neighbours )) = R(k,neighbours);
% end
% 
% 
% 
% 
% 
% 
% %Next add circles! (do not move line back!)
% [all_cycle, description] = f_all_cycle(Adjac);
% short_cycles = cellfun(@(x) length(x), all_cycle);
% 
% for k = 1:numel(all_cycle)
%     for l = 1:length(all_cycle{k})-1
%         B(n-N_bath+k,Index(all_cycle{k}(l), all_cycle{k}(l+1))) = sign(all_cycle{k}(l+1) - all_cycle{k}(l));
%     end
%     B(n-N_bath+k, Index(all_cycle{k}(end), all_cycle{k}(1))) = sign(all_cycle{k}(1) - all_cycle{k}(end));
% end
%        
% pot = zeros(1,size(B,2));
% pot([Index(n-1, 1), Index(1, 4), Index(4,6)]) = [-1, 1, 1];
% B = [B; pot]
% B(:,[5 7 8 9 11 12  13 15]) = []
% B \ [zeros(7,1); 1]






