function mut_info = mutual_information(rho,dimensions)


S = partial_trace_2(rho,dimensions, 2);
B = partial_trace_2(rho,dimensions, 1);

mut_info = entropy(S) + entropy(B) - entropy(rho);