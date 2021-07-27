





bas3 = Basis(3,[],[],[],1)
geom3 = Geometry(3,'f', ones(3))
H_sp3_I = [2,3,2; 3,1,1; 2,1,-1]

H_sp3_N = [2,3,0; 3,1,0; 0,0,-1]

Hub_3_I = Hubbard(geom3, H_sp3_I, [],[],[],bas3,false)
Hub_3_N = Hubbard(geom3, H_sp3_N, [],[],[],bas3,false)

bas2 = Basis(2,[],[],[],1)
geom2 = Geometry(2,'f', ones(2));

H_sp2 = H_sp3_I(1:2,1:2)
Hub_2 = Hubbard(geom2, H_sp2, [],[],[],bas2,false)

full(Hub_3_I.Hamiltonian_mb(ind, ind))
full(Hub_3_N.Hamiltonian_mb(ind, ind))

C = full(Hub_3_I.Hamiltonian_mb(ind,ind))
reshape(sum(reshape(reshape(C(kron(ones(4), eye(2)) == 1), 4,[]).', 2,[]),1),4,[]).'

reshape(sum(reshape(C(kron(eye(4), ones(2)) == 1), 4,[]), 1), 2,[])


% inversion of C = A \otimes eye + eye \otimes B is given by
% (Tr_B (C) - eye(dim_A) * Tr(B) )/dim_B -> A

diff_mat = full(Hub_3_N.Hamiltonian_mb - Hub_3_I.Hamiltonian_mb)
L1 = kron(eye(4), ones(2))
L2 = kron(ones(4), eye(2))