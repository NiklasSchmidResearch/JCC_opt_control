% DP for maximum safe policy
for k=flip(1:N)
    [V_safest(:,k), indices_of_safest_inputs] = max(pagemtimes(T_u, V_safest(:,k+1)), [], 3);
    V_safest(indices_of_unsafe_states, k) = 0;

    % Compute corresponding cost
    T_optimal = zeros(num_vars);
    cost_of_safest_inputs = zeros(num_vars, 1);
    for row_idx = 1:num_vars
        T_optimal(row_idx,:) = T_u(row_idx, :, indices_of_safest_inputs(row_idx));
        cost_of_safest_inputs(row_idx) = stage_cost(row_idx, indices_of_safest_inputs(row_idx));
    end
    C_safest(:,k) = T_optimal*C_safest(:,k+1) + cost_of_safest_inputs;
end