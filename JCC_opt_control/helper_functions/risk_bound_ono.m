mu_k_star_minrisk = zeros(num_vars,N);
for k=flip(1:N)
    [R_bound(:,k), mu_k_star_minrisk(:,k)] = min(pagemtimes(T_u,R_bound(:,k+1)) + indicator_of_A_complement,[],3);
end