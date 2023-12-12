% @authors: Niklas Schmid, Marta Fochesato, Sarah H.Q. Li, Tobias Sutter, John Lygeros @ ETHZ Zurich,
% Automatic Control Laboratory

% Code for paper:
% ### Computing Optimal Joint Chance Constrained Control Policies ###

%% Clear console, workspace, figures
clc
clear all
close all
mycolors = [.8 0.8 0.8; 1 1 1];
% Add subfolders to path
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

%% Settings
disp("Setting parameters.")
% Model parameters
L=40;
r=1;
rng(1234)
M = 10;
sigma = 5;
Delta = 10^(-6);

Mx=round(1.5*L); My=1;  % Discretization of the states
Mu = 5;                 % Discretization of the inputs                                                 
Md = 500;               % Disturbance Samples
Max_Iterations = 30; %100;    % Bound number of iterations

%% Assign states to belong to the safe or unsafe set.
disp("Defining mask.")
% Safe set: states are assigned a value of 2 in the mask 
% Unsafe set: states are assigned a value of 0 in the mask
mask = zeros(Mx,My);

N = 100;                        % Time horizon of the control task
safety_limit_of_states = 13;    % At least 13 units of biomass to be safe
mask(safety_limit_of_states:end) = 2*ones(size(mask(safety_limit_of_states:end)));
current_state = 20;             % Initial state

p_x=0.75;                       % Safety probability

figure(1)
% Plot mask
imagesc(mask)
hold on
scatter(current_state,1);  % Switching x,y because x is rows for us (vertical axis) and y is columns, but matlab plots other way around
title('mask with safe and unsafe states')
%% Compute tranistion kernel
% Define an index to all states for easier reference. Also, keep a list of
% indices belonging to safe and unsafe states. Then, compute transition 
% kernel and stage cost for every state-action pair.
num_vars = Mx*My;                           % Total number of states
stage_cost = zeros(num_vars,Mu);            % State cost assigned every state-action pair
idx_map = zeros(Mx,My);                     % Index belonging to a state at given position
indices_of_safe_states = zeros(Mx,My);      % List of indices of safe states
indices_of_unsafe_states = zeros(Mx,My);    % List of indices of unsafe states

filterMaskAndComputeKernelFish;                 % Compute transition kernel and stage costs

%% Compute the safest and cheapest policy
disp("Computing invariance.")
clims = [0 1];

% Set up value functions for safest and cheapest policy
V_safest = ones(num_vars, N+1); % \bar{V} in paper
V_cheapest = ones(num_vars,N+1);  % \underline{V} in paper
V_safest(indices_of_unsafe_states, N+1) = 0;
V_cheapest(indices_of_unsafe_states, N+1) = 0;
C_safest = zeros(num_vars, N+1);
C_cheapest = zeros(num_vars,N+1);

safest_dp;  % Compute safest policy
cheapest_dp; % Compute cheapest policy


% figure(4)
% V_map = reshape(stage_cost(:,1),Mx,My);
% imagesc(V_map)

% For all three different approaches, at every iteration, store safety and cost
results_of_approaches = zeros(3,Max_Iterations,2); 

%% Approach by Ono et al. (modified due to infeasibility, see paper)
lambda_low = 0;
lambda_up = 10^5;
lambda = lambda_up;

unsafe_indeces_vector = zeros(num_vars,1);
unsafe_indeces_vector(indices_of_unsafe_states) = ones(numel(indices_of_unsafe_states),1);

for iterations=1:Max_Iterations
    % Initialize value function
    J_ono = zeros(size(C_cheapest));
    J_ono(:,N+1)=lambda*unsafe_indeces_vector;

    C_ono = zeros(size(C_cheapest));    % Cost of the policy

    V_ono = zeros(size(V_cheapest));    % Safety of the policy
    V_ono(indices_of_safe_states,N+1) = ones(numel(indices_of_safe_states),1);

    mu_k_star_ono=zeros(num_vars,N);    % Policy
    
    for k=flip(1:N)                     % Run DP recursion
        [J_ono(:,k), mu_k_star_ono(:,k)] = min(pagemtimes(T_u,J_ono(:,k+1)) + reshape(stage_cost,num_vars,1,Mu) + lambda*unsafe_indeces_vector,[],3);
        T_optimal_cost = zeros(num_vars,num_vars);
        for state_idx = 1:num_vars  % Compute cost and safety under the optimal input
            C_ono(state_idx,k) = T_u(state_idx,:,mu_k_star_ono(state_idx,k))*C_ono(:,k+1) + stage_cost(state_idx, mu_k_star_ono(state_idx,k));
            V_ono(state_idx,k) = T_u(state_idx,:,mu_k_star_ono(state_idx,k))*V_ono(:,k+1);
        end
        V_ono(indices_of_unsafe_states, k) = 0;
    end

    % Bisection (deviating from implementation by Ono due to infeasibility)
    if V_ono(current_state,1)<p_x
        lambda_low = lambda
        lambda = 1/2*(lambda_low + lambda_up);
    else
        lambda_up = lambda
        lambda = 1/2*(lambda_low + lambda_up);
        tmp_ono.pol = mu_k_star_ono;
        tmp_ono.C = C_ono(current_state,1);
        tmp_ono.V = V_ono(current_state,1);
    end

    % Store result for current lambda
    results_of_approaches(1,iterations,1) = C_ono(current_state,1);
    results_of_approaches(1,iterations,2) = V_ono(current_state,1);
end

%% Our approach
lambda_low = 0;
lambda_up = (C_safest(current_state,1) - C_cheapest(current_state,1)) / (V_safest(current_state,1) - p_x); %10^5;
lambda = lambda_up; 

V_lambda_up = 0; C_lambda_up = 0;
V_lambda_low = V_cheapest(current_state,1); C_lambda_low = C_cheapest(current_state,1);

not_breaking = true;
iterations=1;
% Outer loop iteration: "Minimize lambda while remaining safe enough"
while not_breaking
    % Initialize value function
    J_our = C_cheapest;
    J_our(indices_of_unsafe_states,:) = J_our(indices_of_unsafe_states,:) + lambda;
    % Value function for cost
    C_our = zeros(size(C_cheapest));
    C_our(indices_of_unsafe_states,:) = C_cheapest(indices_of_unsafe_states,:);
    % Value function for safety
    V_our = zeros(size(V_cheapest));
    V_our(indices_of_safe_states,N+1) = ones(numel(indices_of_safe_states),1);
    mu_k_star_our=zeros(numel(indices_of_safe_states),N);
    % Run DP Recursion
    for k=flip(1:N)
        [J_our(indices_of_safe_states,k), mu_k_star_our(:,k)] = min(pagemtimes(T_u(indices_of_safe_states,:,:),J_our(:,k+1)) + reshape(stage_cost(indices_of_safe_states,:),numel(indices_of_safe_states),1,Mu),[],3);
        for row_idx = 1:numel(indices_of_safe_states) % Compute cost and safety of optimal inputs
            C_our(indices_of_safe_states(row_idx),k) = T_u(indices_of_safe_states(row_idx),:,mu_k_star_our(row_idx,k))*C_our(:,k+1) + stage_cost(indices_of_safe_states(row_idx),mu_k_star_our(row_idx,k));
            V_our(indices_of_safe_states(row_idx),k) = T_u(indices_of_safe_states(row_idx),:,mu_k_star_our(row_idx,k))*V_our(:,k+1);
        end
        V_our(indices_of_unsafe_states, k) = 0;
    end

    % Bisection
    if V_our(current_state,1)<p_x
        lambda_low = lambda
        lambda = 1/2*(lambda_low + lambda_up);
        V_lambda_low = V_our(current_state,1);
        C_lambda_low = C_our(current_state,1);
    else
        lambda_up = lambda
        lambda = 1/2*(lambda_low + lambda_up);
        tmp_our.pol = mu_k_star_our;
        tmp_our.C = C_our(current_state,1);
        tmp_our.V = V_our(current_state,1);
        V_lambda_up = V_our(current_state,1);
        C_lambda_up = C_our(current_state,1);
    end

    % Interpolate to obtain stochastic policy
    p_up = (p_x - V_lambda_low) / (V_lambda_up - V_lambda_low);
    p_low = 1 - p_up;
    C_our_stochastic = C_lambda_low*p_low + C_lambda_up*p_up;
    V_our_stochastic = V_lambda_low*p_low + V_lambda_up*p_up; % Should always equal to p_s

    % Store result for current lambda
    results_of_approaches(2,iterations,1)=C_our_stochastic;
    results_of_approaches(2,iterations,2)=V_our_stochastic;

    suboptimality_bound = p_up*p_low*(lambda_up - lambda_low)*(V_lambda_up - V_lambda_low);
    if suboptimality_bound <= Delta
        fprintf("Breaking after %1.0f iterations because the desired suboptimality is attained.", iterations);
        not_breaking = false;
    end
    iterations = iterations + 1;
end

%% Plot safety cost pairs obtained over the outer loop iterations
% figure(20)
% plot(results_of_approaches(1,:,1),results_of_approaches(1,:,2))
% hold on
% plot(results_of_approaches(2,:,1),results_of_approaches(2,:,2))
% legend('Ono et al.','Ours')
% xlabel('Cost')
% ylabel('Safety')
fprintf('Ono et al, Our approach (cost,safety), (%4.2f,%4.2f),(%4.2f,%4.2f)',tmp_ono.C,tmp_ono.V,tmp_our.C,tmp_our.V)

% For simplicity we just stick to the policy associated to lambda_up
% instead of the mixed policy. 
mu_k_star_ono = tmp_ono.pol;
mu_k_star_our = tmp_our.pol;

%% Run simulation
figure(14)
hold on

mu_k_star_our_ext = zeros(size(mu_k_star));
mu_k_star_our_ext(indices_of_safe_states,:) = mu_k_star_our;

success_rate = 0;
total_attempts=5000;        % Number of Monte-Carlo runs
cost_avg = 0;

for attempt=1:total_attempts
    cost_of_attempt = 0;
    success = 1;
    x_idx = zeros(N+1,1);
    x_idx(1) = current_state;
    x = zeros(N+1,1);
    x(1) = current_state -1;

    for k=1:N
        if success==1
            u = mu_k_star_our_ext(x_idx(k),k);
        else
            u = mu_k_star(x_idx(k),k);
        end

        %u = mu_k_star_ono(x_idx(k),k);   % Uncomment to evaluate Ono's policy

        %u = mu_k_star(x_idx(k),k);       % Uncomment to evaluate cheapest policy
        
        % Run simulation step
        d = (u-1) / Mu;
        v = normrnd(0.2, 0.1^2);
        v = boundValue(v,[0,1]);
        gamma = normrnd(1, 0.6^2);
        delta = normrnd(1.1, 0.2^2);
        delta = boundValue(delta,[0,inf]);
        x_k = x(k);
        C = delta*d*M*x_k/L;
        C = boundValue(C,[0,delta*d*M]);
        K = 1 ./ (1 + exp(-(x_k-20)/sigma) );
        R = r*x_k*(1- x_k/L);
        R = max(R,0);
        x(k+1) = (1-v)*x_k + gamma*R*K - C;
        x(k+1) = round(max(x(k+1),0));
        x_idx(k+1) = x(k+1)+1;
    
        % Bound values of state 
        if x_idx(k+1)<=0
            x_idx(k+1)=1;
        end
        if x_idx(k+1)>Mx
            x_idx(k+1) = Mx;
        end

        % Check safety
        x(k+1) = x_idx(k+1) - 1;
        if abs(mask(x_idx(k+1)))<1
            success = 0;
        end

        cost_of_attempt = cost_of_attempt + stage_cost(x_idx(k),u);
    end
    cost_avg = cost_avg + cost_of_attempt/total_attempts;
    success_rate = success_rate + success/total_attempts;
    plot(x,'LineWidth',1)
end
xlabel('Time')
ylabel('Units of Biomass')
title('Monte-Carlo Simulation')
cost_avg
safety =success_rate;
safety

%% Surface plots of policies
fontsize = 20;
colorrange = linspace(0,1,1000)';
colorrangeZero =zeros(size(colorrange));
colorrangeOne = ones(size(colorrange))*1;
%mycolors = [0 0 0.8; 0 0.8 0; 0.8 0.8 0];
mycolors = [colorrangeZero flip(colorrange) colorrange ; colorrange colorrangeZero flip(colorrange)];
mycolors = flip([colorrangeOne colorrange colorrange; flip(colorrange) flip(colorrange) colorrangeOne ]);
mycolors = flip([colorrangeOne*0, flip(colorrange)*1 + 0, colorrangeOne*0 ]);
mycolors = flip([flip(colorrange)*1, flip(colorrange)*.5 + 0.5, flip(colorrange)*1 ]);
%mycolors = [colorrange colorrange flip(colorrange); colorrangeOne flip(colorrange) colorrangeZero];
%mycolors = [flip(colorrange) flip(colorrange) flip(colorrange)];
figure(33)
a = linspace(safety_limit_of_states,N-1,size(mu_k_star_ono(indices_of_safe_states,:),2));
b = linspace(safety_limit_of_states,N-1,size(mu_k_star_ono(indices_of_safe_states,:),1));
imagesc(a,b,mu_k_star_ono(indices_of_safe_states,:)/Mu)
colormap(mycolors);       
caxis([-1, 1]);    
ax = gca;
set(gca, 'FontName', 'times')
ax.FontSize = fontsize; 
ylabel('State', 'FontSize', fontsize)
xlabel('Timestep k', 'FontSize', fontsize)
set(gca,'FontName','Times Roman','Fontsize',fontsize)
title('Policy k=0 by Ono et al.')
figure(34)
a = linspace(safety_limit_of_states,N-1,size(mu_k_star_our_ext(indices_of_safe_states,:),2));
b = linspace(safety_limit_of_states,N-1,size(mu_k_star_our_ext(indices_of_safe_states,:),1));
imagesc(a,b,mu_k_star_our_ext(indices_of_safe_states,:)/Mu)
colormap(mycolors);      
caxis([-1, 1]);   
ax = gca;
set(gca, 'FontName', 'times')
ax.FontSize = fontsize; 
ylabel('State', 'FontSize', fontsize)
xlabel('Timestep k', 'FontSize', fontsize)   
set(gca,'FontName','Times Roman','Fontsize',fontsize)
title('Policy k=0 by Our Approach, b_k=1')
figure(35)
diff = mu_k_star_our_ext(indices_of_safe_states,:) - mu_k_star_ono(indices_of_safe_states,:);
a = linspace(safety_limit_of_states,N-1,size(diff,2));
b = linspace(safety_limit_of_states,N-1,size(diff,1));
imagesc(a,b,diff/Mu)
colormap(mycolors);    
caxis([-1, 1]);         
colorbar
ax = gca;
ax.FontSize = fontsize; 
ylabel('State', 'FontSize', fontsize)
xlabel('Timestep k', 'FontSize', fontsize)
set(gca, 'FontName', 'times')
set(gca,'FontName','Times Roman','Fontsize',fontsize)
title('Difference of Policies')

