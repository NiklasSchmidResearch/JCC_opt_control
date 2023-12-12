

% Settings of simulation
disp("Setting parameters.")
%discrete_states_per_dimension = [50,50]; % 50 x 50 discretization in the x- and y-direction
%discrete_actions_per_dimension = [8]; % 8 inputs as discretization of U
Mx=50; My=50;                                                                               % REMOVE 
Mu = 8; % Discretization of the input set                                                   % REMOVE
Max_Iterations_Ono = 30;  % For the approach in Ono we set some high number on the maximum iterations
scenario = 2;       % Simulated scenario (see quadcopter_example.m)
Delta = 10e-6; % Maximum suboptimality

COMPUTE_PARETO_FRONT = true;


% Some settings for policy plots
mycolors = [.8 0.8 0.8; 1 1 1];
quivercolors = [0.7 1 0.7; 1 1 1];
quiver_spacing=5;
quiver_arrowscaling = 0.4;
headWidth =3;  % 1/10 of annotation
headLength=1.5;  % 1/10 of annotation
LineLength = 1.3; % same as annotation
linestyle = '-k2';

% Set random seed
rng(1234)
