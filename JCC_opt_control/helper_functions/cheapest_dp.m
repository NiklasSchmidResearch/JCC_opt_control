mu_k_star = zeros(num_vars,N);
for k=flip(1:N)
    [C_cheapest(:,k), mu_k_star(:,k)] = min(pagemtimes(T_u,C_cheapest(:,k+1)) + reshape(stage_cost,num_vars,1,Mu),[],3);
    T_optimal = zeros(num_vars);
    for row_idx = 1:num_vars
        T_optimal(row_idx,:) = T_u(row_idx,:,mu_k_star(row_idx,k));
    end
    V_cheapest(:,k) = T_optimal*V_cheapest(:,k+1);
    V_cheapest(indices_of_unsafe_states, k) = 0;
end

    
