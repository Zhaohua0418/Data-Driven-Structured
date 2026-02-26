clear
rng(45, 'twister');
nx = 3;
nu = 2;
ny = 5;
nd = 3;

S = [1 1 0; 
    0 1 1];

% epsilon list
epsilon_list = [0.05; 0.1; 0.15; 0.2; 0.5];

% T = 6;
T = 10;
% T = 20;
load("H2_traj_T_20.mat");

% system parameters
sys = struct( ...
    'A', [ -0.4095     0.4036   -0.0874
            0.5154    -0.0815    0.1069
            1.6715     0.7718   -0.3376], ...
    'B', [  0          0
           -0.6359    -0.1098
           -0.0325     2.2795],...
    'C', [ eye(nx);...
           zeros(nu, nx)],...
    'D', [ zeros(nx,nu);...
           eye(nu)],...
    'G', eye(nx),...
    'nx',nx,'nu',nu,'ny',ny,'nd',nd);

%%
h2_module = H2_module_my_structured(sys, S);

%% model-based and data-driven control
sol_table = cell(1, 1+length(epsilon_list));
h2_table = NaN*zeros(1, 1+length(epsilon_list));
sol_table{1} = h2_module.model_based_H2(); % cell一定要用{}
h2_table(1) = sol_table{1}.h2;
for i=1:length(epsilon_list)
    try
        sol_table{i+1} = h2_module.data_driven_H2(traj{i});
        h2_table(i+1) = sol_table{i+1}.h2;
    catch % catch infeasibility
        sol_table{i} = NaN;
    end
end