clear all
rng(45, 'twister');
cvx_solver mosek
nx = 3;
nu = 2;
ny = 5;
nd = 3;
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

S = [1 1 0; 
    0 1 1];

% epsilon fixed to 0.1


fname = sprintf('H2_info_struct_epsilon_0.1_T_200.mat');
load(fname);

H2_info_struct_temp = H2_info_struct;
%%
h2_my_module = H2_module_my_structured(sys, S);
h2_jared_module = H2_module_jared_structured(sys, S);
my_running_time = [];
jared_running_time = [];

%% data-driven control
for t = 5:5:200
    H2_info_struct_temp.T = t;
    H2_info_struct_temp.X = H2_info_struct.X(:,1:t+1);
    H2_info_struct_temp.U = H2_info_struct.U(:,1:t);
    start_time = tic;
    try
        h2_my_module.data_driven_H2(H2_info_struct_temp);
    catch
    end
    elapsed_time = toc(start_time);
    my_running_time = [my_running_time elapsed_time];
    start_time = tic;
    try
    h2_jared_module.data_driven_H2(H2_info_struct_temp);
    catch
    end
    elapsed_time = toc(start_time);
    jared_running_time = [jared_running_time elapsed_time];
end
