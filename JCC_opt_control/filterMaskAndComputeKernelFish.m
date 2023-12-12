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

tic
for i_x=1:Mx   % For all x positions
    i_y = 1;
    idx_i = transformMatrixIdxToVectorIdx([Mx,My],[i_x,i_y]);
    for u=1:Mu   % try all inpus
        for sampling_iteration=1:Md
            % System model
            d = (u-1) / Mu;
    
            v = normrnd(0.2, 0.1^2);
            v = boundValue(v,[0,1]);
            gamma = normrnd(1, 0.6^2);
            delta = normrnd(1.1, 0.2^2);
            delta = boundValue(delta,[0,inf]);
        
            x_k = i_x - 1;
            % Fishing:
            C = delta*d*M*x_k/L;
            C = boundValue(C,[0,delta*d*M]);
    
            K = 1 ./ (1 + exp(-(x_k-20)/sigma) );
            R = r*x_k*(1- x_k/L);
            R = max(R,0);
            x_kp = (1-v)*x_k + gamma*R*K - C;
            x_kp = round(max(x_kp,0));
    
            j_x = x_kp + 1;
            j_y = 1;
    
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
    
            idx_j = transformMatrixIdxToVectorIdx([Mx,My],[j_x,j_y]);
            % Compute indices, increase by Md (add loop), stage cost is
            % negative catch
    
            % Put transition probabilities in transition matrix
            T_u(idx_i,idx_j,u) = T_u(idx_i,idx_j,u) + 1/Md;
    
            % Stage cost is quadratic. Zero in the middle.
            stage_cost(idx_i,u) = - C;
            terminal_cost(idx_i) = 0;
        end
    end
end
toc

