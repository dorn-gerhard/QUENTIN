
n = 35;
A = rand(150,240);
B = eye(n);

% ----------------- eye_kron test -----------------
tic
for k = 1:100
    C = kron(B,A);
end
toc
% 11 sec


profile on
tic
for k = 1:100
    D = eye_kron(n,A);
end
toc
profile viewer

% 0.1 sec

% ----------------- kron_eye test -----------------
tic
for k = 1:100
    C = kron(A,B);
end
toc
% 11 sec


profile on
tic
for k = 1:100
    D = kron_eye(A,n);
end
toc
% 0.9 sec
profile viewer


% ----------------- combindes kron_ones test -----------------
A = rand(80,100);

B = ones(size(A));




% ----------------- ones_kron test -----------------


tic
for k = 1:100
    C = kron(conj(A), ones(size(A))) + kron(ones(size(A)), A);
    disp(k)
end
toc
% 48 sec


tic
profile on
for k = 1:100
    D = ones_kron(A);
    disp(k)
end
profile viewer
toc
% 40 sec



