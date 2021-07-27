function S = random_rho(n)


S = rand(n)-0.5 + (rand(n)-0.5)*1i;
S = (S+S')/2;
[ev,ew] = eig(S);
S = ev*abs(ew) * ev';
S = S./trace(S);