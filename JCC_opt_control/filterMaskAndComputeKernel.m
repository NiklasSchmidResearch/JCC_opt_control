%% Iterate through all values in the safe set and assign indeces to every representative point in the safe set

% Assign a number (index) to each state. Store lists of indeces belonging
% to unsafe or target sets.
%idx_map = reshape(1:num_vars, Mx, My);
indices_of_safe_states = find(abs(mask-2)<0.01);
indices_of_unsafe_states = find(abs(mask)<0.01);

%%% Run an exhaustive simulation to compute the transition kernel
disp("Building transition kernel by exhaustive simulation. This may take some time.")
T_u = zeros(num_vars,num_vars,Mu);      % Transition kernel
c_u = zeros(num_vars,Mu);               % Stage cost function

% Start by defining grid points over a space twice as large as the state
% space
[grid_x,grid_y] = meshgrid(linspace(-Mx,Mx,2*Mx+1),linspace(-My,My,2*My+1));
grid_over_space = reshape([grid_y,grid_x],[],2)'; 
% grid_over_space goes over x in the inner, y in the outer loop, in the end
% does not matter because our disturbance distribution is symmetric

% Compute disturbance distribution once, since we make it state-invariant
normal_distribution_true = zeros(2*Mx+1,2*My+1,1);       % Distribution of noise
for k=1:size(grid_over_space,2) % Go through all grid points and compute probability
    matrix_idx = transformVectorIdxToMatrixIdx(size(normal_distribution_true),k)';
    normal_distribution_true(matrix_idx(1),matrix_idx(2)) = mvnpdf(grid_over_space(:,k),[0;0],Sigma_disturbance);
end

Md = 100000;
normal_distribution = zeros(2*Mx+1,2*My+1,1);       % Distribution of noise
for sampling_idx=1:Md
    d = mvnrnd([0;0],Sigma_disturbance);
    d1 = round(d(1));
    d2 = round(d(2));
    normal_distribution(d1 + 1 + Mx, d2 + 1 + My) = normal_distribution(d1 + 1 + Mx, d2 + 1 + My) + 1/Md;
end



% Compute maximum and plot distribution of disturbance
figure(123)
imagesc(normal_distribution)
title('Noise Distribution')

tic
for i_x=1:Mx   % For all x positions
    for i_y=1:My % For all y positions
        idx_i = transformMatrixIdxToVectorIdx([Mx,My],[i_x,i_y]);
        for u=1:Mu   % try all inpus
            % Stage cost is quadratic. Zero in the middle.
            stage_cost(idx_i,u) = norm([i_x - Mx/2 - 0.5; i_y - My/2 - 0.5],2)^2;
            terminal_cost(idx_i) = norm([i_x - Mx/2 - 0.5; i_y - My/2 - 0.5],2)^2;

            % System model
            psi = u / Mu*2*pi; 
            j_x = i_x + round(3*Mx/50*cos(psi));
            j_y = i_y + round(3*My/50*sin(psi));
            

            % We are bounded to end up in the state space
            if j_x<=0
                j_x=1;
            end
            if j_x>Mx
                j_x = Mx;
            end
            if j_y<=0
                j_y=1;
            end
            if j_y>My
                j_y=My;
            end              

            % Based on disturbance distribution, obtain stochastic model
            % simply by placing the distribution the state we ended up in
            transition_distribution = normal_distribution(1 + Mx + 1 - j_x: 2*Mx + 1 - j_x, ...
                1 + My + 1 - j_y : 2*My + 1 - j_y);
            % Normalize distribution
            transition_distribution = transition_distribution/sum(sum(transition_distribution));
            % Put transition probabilities in transition matrix
            T_u(idx_i,:,u) = reshape(transition_distribution, num_vars, 1);
        end
    end
end
toc

