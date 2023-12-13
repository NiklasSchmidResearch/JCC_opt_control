% @authors: Niklas Schmid, Marta Fochesato, Sarah H.Q. Li, Tobias Sutter, John Lygeros @ ETHZ Zurich,
% Automatic Control Laboratory

% Code for paper:
% ### Computing Optimal Joint Chance Constrained Control Policies ###

% We use the arrow3 script by Tom Davis from the matlab add-ons to plot the
% policies

%% Clear console, workspace, figures
clc
clear all
close all

% Add subfolders to path
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

%% Load setttings
quadcopter_example_settings;

%% Assign states to belong to the safe or unsafe set.
% Safe set: states are assigned a value of 2 in the mask 
% Unsafe set: states are assigned a value of 0 in the mask

%W.l.o.g. we shifted the state space by (+25,+25) to range from [0,50]x[0,50]
disp("Defining mask.")
mask = zeros(Mx,My);
if scenario == 1
    N = 20; % Time horizon of the control task
    Sigma_disturbance = 5*eye(2);                                   % Noise variance
    range_inner_circle=5/16;
    mask(round(Mx*1/8):Mx - round(Mx*1/8),round(My*1/8):My - round(My*1/8)) = 2*ones(size(mask(round(Mx*1/8):Mx - round(Mx*1/8),round(My*1/8):My - round(My*1/8))));
    mask(round(Mx*range_inner_circle):Mx - round(Mx*range_inner_circle),round(My*range_inner_circle):My - round(My*range_inner_circle)) = zeros(Mx - 2*round(Mx*range_inner_circle) +1 , My - 2*round(My*range_inner_circle) +1);
    x_0_idx=transformMatrixIdxToVectorIdx([Mx,My],[round(Mx*6/8),round(My*6/8)]);
    p_s = 0.6;
elseif scenario == 2
    N = 20; % Time horizon of the control task
    Sigma_disturbance = 5*eye(2);                                   % Noise variance
    mask = 2*ones(Mx,My);
    mask(round(Mx*4.5/8):Mx - round(Mx*1.5/8),round(My*4/8):My - round(My*1/8)) = 0*mask(round(Mx*4.5/8):Mx - round(Mx*1.5/8),round(My*4/8):My - round(My*1/8));
    x_0_idx=transformMatrixIdxToVectorIdx([Mx,My],[round(Mx*7/8),round(My*7/8)]);
    p_s = 0.9;
end

% Plot mask and initial state
figure(1)
imagesc(mask)                               % Mask
title('Mask')
hold on
x0 = transformVectorIdxToMatrixIdx([Mx,My],x_0_idx);
scatter(x0(1,1),x0(2,1),90,'white','filled','o','LineWidth',3,'MarkerEdgeColor','black','SizeData',400)
xlabel('x-coordinate')
ylabel('y-coordinate')
set(gca, 'XTick', [1:10:50,50], 'XTickLabel', [-25:10:25,25]) % 10 ticks 
set(gca, 'YTick', [1:10:50,50], 'YTickLabel', [-25:10:25,25]) % 10 ticks 
set(gca,'FontName','Times Roman','Fontsize',25)
colormap(mycolors);
set(gca,'FontName','Times Roman','Fontsize',25)

%% Compute tranistion kernel
% Define an index for all states for easier reference. Also, keep a list of
% indices belonging to safe and unsafe states. Then, compute transition 
% kernel and stage cost for every state-action pair. In the implementation,
% we treat the state space as one-dimensional by stacking all states in a
% vector.
num_vars = Mx*My;                           % Total number of states
stage_cost = zeros(num_vars,Mu);            % State cost assigned every state-action pair
idx_map = zeros(Mx,My);                     % Index belonging to a state at given position
indices_of_safe_states = zeros(Mx,My);      % List of indices of safe states
indices_of_unsafe_states = zeros(Mx,My);    % List of indices of unsafe states

filterMaskAndComputeKernel;                 % Compute transition kernel and stage costs

%% Compute the safest and cheapest policy
disp("Checking triviality and feasibility.")
clims = [0 1];

% Set up value functions for safest (denoted with overline in paper) and cheapest (denoted with underline in paper) policy
% Safety is set identity for last time-step
V_safest = ones(num_vars, N+1); 
V_cheapest = ones(num_vars,N+1);  
V_safest(indices_of_unsafe_states, N+1) = 0;
V_cheapest(indices_of_unsafe_states, N+1) = 0;
% Cost is zero at last time-step since we do not assign a terminal cost for
% the quadcopter example
C_safest = zeros(num_vars, N+1);
C_cheapest = zeros(num_vars,N+1);

% Run DP recursion for safest and cheapest policies
safest_dp;  % Compute safest policy
cheapest_dp; % Compute cheapest policy

% Plot the stage cost over the state space for some input
figure(4)
V_map = reshape(stage_cost(:,1),Mx,My);
imagesc(V_map)
title('cost function')

% We will test three CCOC approaches in the following and store the
% performance of the resulting policies in the following variable.
results_of_approaches = zeros(3,Max_Iterations_Ono,2); 

%% Compute feasibility via Ono's approach and via our approach
% Create a vector with ones at unsafe states
indicator_of_A_complement = zeros(num_vars,1);
indicator_of_A_complement(indices_of_unsafe_states) = ones(numel(indices_of_unsafe_states),1);
R_bound = zeros(size(V_cheapest));
risk_bound_ono;
fprintf(['The bound in Ono returns a minimum risk of %1.4f, meaning a safety' ...
    ' of %1.4f, while we require %1.4f and could achieve %1.4f\n'], ...
    R_bound(x_0_idx,1),1-R_bound(x_0_idx,1),p_s,V_safest(x_0_idx,1));

%% Approach by Ono et al.
disp("Running the approach by Ono.")
% Initialize Bisection algorithm
lambda_low = 0;
lambda_up = 10^5;   % Choose some high value here
lambda = lambda_up;


% Outer loop iteration: "Minimize lambda while remaining safe enough"
for iterations=1:Max_Iterations_Ono
    % J: Value function of DP recursion, C, V: Cost, safety of resulting,
    % R: Risk according to bound in Ono's paper
    J_ono = zeros(size(C_cheapest));
    J_ono(:,N+1)=lambda*indicator_of_A_complement;

    C_ono = zeros(size(C_cheapest));
    
    V_ono = zeros(size(V_cheapest));
    V_ono(indices_of_safe_states,N+1) = ones(numel(indices_of_safe_states),1);

    R_ono = zeros(size(V_cheapest));
    R_ono(indices_of_unsafe_states,N+1) = ones(numel(indices_of_unsafe_states),1);  % In ono variant set unsafe states to one
    
    % Store the optimal policy under the variable
    pi_star_ono=zeros(num_vars,N);    
    
    % Inner loop DP recursion
    for k=flip(1:N)
        % J = min_{mu_k} T_{mu_k}*J(x_k) + l(x,u) + lambda*1_{A^c}, where
        % T_{mu_k} is the transition kernel under mapping mu_k.
        [J_ono(:,k), pi_star_ono(:,k)] = min(pagemtimes(T_u,J_ono(:,k+1)) + reshape(stage_cost,num_vars,1,Mu) + lambda*indicator_of_A_complement,[],3);
        % Compute the cost, safety and risk associated with this policy
        T_optimal_cost = zeros(num_vars,num_vars);
        for state_idx = 1:num_vars
            C_ono(state_idx,k) = T_u(state_idx,:,pi_star_ono(state_idx,k))*C_ono(:,k+1) + stage_cost(state_idx, pi_star_ono(state_idx,k));
            V_ono(state_idx,k) = T_u(state_idx,:,pi_star_ono(state_idx,k))*V_ono(:,k+1);
            R_ono(state_idx,k) = T_u(state_idx,:,pi_star_ono(state_idx,k))*R_ono(:,k+1) + indicator_of_A_complement(state_idx); % in ono approach keep adding 1 if outside
        end
        V_ono(indices_of_unsafe_states, k) = 0; % Apply indicator function for safety
    end

    % Do bisection. In Ono this would have been done using the risk R, not
    % the true safety V. We instead compute the actual safety here to
    % restrict our comparison to the DP recursions.
    %if R_ono(x_0_idx,1)>1-p_s % <- Use this conidition to be in line with Ono 
    if V_ono(x_0_idx,1)<p_s
        lambda_low = lambda
        lambda = 1/2*(lambda_low + lambda_up);
    else
        lambda_up = lambda
        lambda = 1/2*(lambda_low + lambda_up);
    end

    % Store result for current lambda
    results_of_approaches(1,iterations,1) = C_ono(x_0_idx,1);
    results_of_approaches(1,iterations,2) = V_ono(x_0_idx,1);
    results_of_approaches(3,iterations,1) = C_ono(x_0_idx,1);
    results_of_approaches(3,iterations,2) = R_ono(x_0_idx,1);
end

%% Our approach
disp("Running our approach.")
% Initialize Bisection algorithm
lambda_low = 0;
lambda_up = (C_safest(x_0_idx,1) - C_cheapest(x_0_idx,1)) / (V_safest(x_0_idx,1) - p_s); %10^5;
lambda = lambda_up;

V_lambda_up = 0; C_lambda_up = 0;
V_lambda_low = V_cheapest(x_0_idx,1); C_lambda_low = C_cheapest(x_0_idx,1);



mu_k_star_our_up=zeros(numel(indices_of_safe_states),N);
mu_k_star_our_low=zeros(numel(indices_of_safe_states),N);

iterations=1;
not_breaking = true;
% Outer loop iteration: "Minimize lambda while remaining safe enough"
while not_breaking
    % J: Value function of DP recursion, C, V: Cost, safety of resulting,
    J_our = C_cheapest;
    J_our(indices_of_safe_states,:) = J_our(indices_of_safe_states,:) - lambda;
    C_our = zeros(size(C_cheapest));
    C_our(indices_of_unsafe_states,:) = C_cheapest(indices_of_unsafe_states,:);
    V_our = zeros(size(V_cheapest));
    V_our(indices_of_safe_states,N+1) = ones(numel(indices_of_safe_states),1);
    % Store the optimal policy under the variable
    mu_k_star_our=zeros(numel(indices_of_safe_states),N);
    % Inner loop DP recursion
    for k=flip(1:N)
        [J_our(indices_of_safe_states,k), mu_k_star_our(:,k)] = min(pagemtimes(T_u(indices_of_safe_states,:,:),J_our(:,k+1)) + reshape(stage_cost(indices_of_safe_states,:),numel(indices_of_safe_states),1,Mu),[],3);
        for row_idx = 1:numel(indices_of_safe_states)
            % Compute the cost, safety and risk associated with this policy
            C_our(indices_of_safe_states(row_idx),k) = T_u(indices_of_safe_states(row_idx),:,mu_k_star_our(row_idx,k))*C_our(:,k+1) + stage_cost(indices_of_safe_states(row_idx),mu_k_star_our(row_idx,k));
            V_our(indices_of_safe_states(row_idx),k) = T_u(indices_of_safe_states(row_idx),:,mu_k_star_our(row_idx,k))*V_our(:,k+1);
        end
        V_our(indices_of_unsafe_states, k) = 0;
    end
    
  % Do bisection. 
    if V_our(x_0_idx,1)<p_s
        lambda_low = lambda
        lambda = 1/2*(lambda_low + lambda_up);
        V_lambda_low = V_our(x_0_idx,1);
        C_lambda_low = C_our(x_0_idx,1);
        mu_k_star_our_low = mu_k_star_our;
    else
        lambda_up = lambda
        lambda = 1/2*(lambda_low + lambda_up);
        V_lambda_up = V_our(x_0_idx,1);
        C_lambda_up = C_our(x_0_idx,1);
        mu_k_star_our_up = mu_k_star_our;
    end

    % Interpolate to obtain stochastic policy
    p_up = (p_s - V_lambda_low) / (V_lambda_up - V_lambda_low);
    p_low = 1 - p_up;
    C_our_stochastic = C_lambda_low*p_low + C_lambda_up*p_up;
    V_our_stochastic = V_lambda_low*p_low + V_lambda_up*p_up; % Should always equal to p_s

    % Store result for current lambda
    results_of_approaches(2,iterations,1)=C_our_stochastic;
    results_of_approaches(2,iterations,2)=V_our_stochastic;

%     V_lambda_low
%     V_lambda_up
%     p_up*p_low
    suboptimality_bound = p_up*p_low*(lambda_up - lambda_low)*(V_lambda_up - V_lambda_low);
    if suboptimality_bound <= Delta
        fprintf("Breaking after %1.0f iterations because the desired suboptimality is attained.", iterations);
        not_breaking = false;
    end

    if iterations==1
        p_up = (p_s - V_cheapest(x_0_idx,1)) / (V_lambda_up - V_cheapest(x_0_idx,1));
        p_low = 1-p_up;
        K = ceil(   	-log2( Delta / (0.25*(lambda_up - lambda_low)) )   );
        fprintf("This will take at maximum %1.0f iterations.", K);
    end

  
    if iterations > K 
        not_breaking = false;
    end
    iterations = iterations + 1;
end


%% Plot the policies and performance of both approaches
% Plot policies
quadcopter_plot_policies;

% Plot performances
figure(20)
plot(results_of_approaches(1,:,1),results_of_approaches(1,:,2))
hold on
plot(results_of_approaches(2,:,1),results_of_approaches(2,:,2))
plot(results_of_approaches(1,:,1),results_of_approaches(1,:,2))
legend('Ono Approach','Our Approach')
title('Pairs of safety and cost obtained during the recursions')

fprintf("Performance evaluation. Ono (C,1-R,V): (%1.4f,%1.4f,%1.4f). Our (C,V): (%1.4f,%1.4f)\n",C_ono(x_0_idx,1),1-R_ono(x_0_idx,1),V_ono(x_0_idx,1),C_our(x_0_idx,1),V_our(x_0_idx,1))


%% Simulate using policies
quadcopter_simulate_policies;

%% Plot Pareto fronts by swiping through lambda
if COMPUTE_PARETO_FRONT
    disp("Computing Pareto fronts. This may take a while. Starting with Ono.\n")
    Iterations_over_lambda = 100;
    lambda_list = logspace(2,6,Iterations_over_lambda);
    
    % Ono approach
    for iterations=1:Iterations_over_lambda
        lambda = lambda_list(iterations);
    
        J_ono = zeros(size(C_cheapest));
        J_ono(:,N+1)=lambda*indicator_of_A_complement;
        C_ono = zeros(size(C_cheapest));
        V_ono = zeros(size(V_cheapest));
        V_ono(indices_of_safe_states,N+1) = ones(numel(indices_of_safe_states),1);
        pi_star_ono=zeros(num_vars,N);
        
        for k=flip(1:N)
            [J_ono(:,k), pi_star_ono(:,k)] = min(pagemtimes(T_u,J_ono(:,k+1)) + reshape(stage_cost,num_vars,1,Mu) + lambda*indicator_of_A_complement,[],3);
            T_optimal_cost = zeros(num_vars,num_vars);
            for state_idx = 1:num_vars
                C_ono(state_idx,k) = T_u(state_idx,:,pi_star_ono(state_idx,k))*C_ono(:,k+1) + stage_cost(state_idx, pi_star_ono(state_idx,k));
                V_ono(state_idx,k) = T_u(state_idx,:,pi_star_ono(state_idx,k))*V_ono(:,k+1);
                R_ono(state_idx,k) = T_u(state_idx,:,pi_star_ono(state_idx,k))*R_ono(:,k+1) + indicator_of_A_complement(state_idx); % in ono approach keep adding 1 if outside
            end
            V_ono(indices_of_unsafe_states, k) = 0;
        end
    
        % Store result for current lambda
        gridded_results(1,iterations,1) = C_ono(x_0_idx,1);
        gridded_results(1,iterations,2) = V_ono(x_0_idx,1);
        gridded_results(3,iterations,1) = C_ono(x_0_idx,1);
        gridded_results(3,iterations,2) = R_ono(x_0_idx,1);
    end
    disp("Computing Pareto front for our approach.\n")
    % Our approach
    for iterations=1:Iterations_over_lambda
        lambda = lambda_list(iterations);
        J_our = C_cheapest;
        J_our(indices_of_safe_states,:) = J_our(indices_of_safe_states,:) - lambda;
        C_our = zeros(size(C_cheapest));
        C_our(indices_of_unsafe_states,:) = C_cheapest(indices_of_unsafe_states,:);
        V_our = zeros(size(V_cheapest));
        V_our(indices_of_safe_states,N+1) = ones(numel(indices_of_safe_states),1);
        mu_k_star_our=zeros(numel(indices_of_safe_states),N);
        for k=flip(1:N)
            [J_our(indices_of_safe_states,k), mu_k_star_our(:,k)] = min(pagemtimes(T_u(indices_of_safe_states,:,:),J_our(:,k+1)) + reshape(stage_cost(indices_of_safe_states,:),numel(indices_of_safe_states),1,Mu),[],3);
            for row_idx = 1:numel(indices_of_safe_states)
                C_our(indices_of_safe_states(row_idx),k) = T_u(indices_of_safe_states(row_idx),:,mu_k_star_our(row_idx,k))*C_our(:,k+1) + stage_cost(indices_of_safe_states(row_idx),mu_k_star_our(row_idx,k));
                V_our(indices_of_safe_states(row_idx),k) = T_u(indices_of_safe_states(row_idx),:,mu_k_star_our(row_idx,k))*V_our(:,k+1);
            end
            V_our(indices_of_unsafe_states, k) = 0;
        end
    
        gridded_results(2,iterations,1)=C_our(x_0_idx,1);
        gridded_results(2,iterations,2)=V_our(x_0_idx,1);
    end
    
    %% Plot the pareto fronts
    figure(220)
    plot(gridded_results(1,:,2),gridded_results(1,:,1))
    hold on
    plot(gridded_results(2,:,2),gridded_results(2,:,1))
    hold on
    indices = (1-gridded_results(3,:,2)>=-0.1);
    plot(1-gridded_results(3,indices,2),gridded_results(3,indices,1))
    legend('Ono''s approach','Our approach','Ono''s approach modified')
    title('Pareto fronts, safety vs cost')
    
    % and print them into the console for the latex plot
    for i = 1:Iterations_over_lambda
        fprintf("(%1.4f,%1.4f)",gridded_results(1,i,2),gridded_results(1,i,1))
    end
    disp("\n")
    for i = 1:Iterations_over_lambda
        fprintf("(%1.4f,%1.4f)",gridded_results(2,i,2),gridded_results(2,i,1))
    end
    disp("\n")
    for i = 1:numel(indices)
        if indices(i)
            fprintf("(%1.4f,%1.4f)",1-gridded_results(3,i,2),gridded_results(3,i,1))
        end
    end
end
