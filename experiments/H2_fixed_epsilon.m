clear all
rng(45, 'twister');
cvx_solver mosek
nx = 3;
nu = 2;
ny = 5;
nd = 3;

S = [1 1 0; 
    0 1 1];

% epsilon fixed to 0.1


T_list = [6; 10; 15];
t = 15;

fname = sprintf('H2_info_struct_epsilon_0.1_T_20.mat');
load(fname);
% to be consistent with jared's code

H2_info_struct_copy = H2_info_struct;

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


%% data-driven control
keywords = ["jared_unstructured_energy_bound","my_unstructured_instantaneous_bound", "jared_structured_energy_bound", "my_structured_energy_bound","jared_structured_instantaneous_bound","my_structured_instantaneous_bound"];
classes = {H2_module_jared_unstructured(sys),H2_module_my_unstructured(sys),H2_module_jared_structured(sys, S),H2_module_my_structured_energy_bound(sys, S),H2_module_jared_structured_instantaneous_bound(sys, S),H2_module_my_structured(sys, S)};
sol_table = cell(length(keywords), 1+length(T_list));
h2_table = NaN*zeros(length(keywords),  1+length(T_list));
for i=1:length(keywords)
    sol_table{i,1} = classes{i}.model_based_H2(); % cell一定要用{}
    h2_table(i,1) = sol_table{i,1}.h2;
    for j=1:length(T_list)
        H2_info_struct_copy.T = T_list(j);
        H2_info_struct_copy.X = H2_info_struct.X(:,1:T_list(j)+1);
        H2_info_struct_copy.U = H2_info_struct.U(:,1:T_list(j));
        try
            sol_table{i,j+1} = classes{i}.data_driven_H2(H2_info_struct_copy);
            h2_table(i,j+1) = sol_table{i,j+1}.h2;
        catch % catch infeasibility
            sol_table{i,j+1} = NaN;
        end
    end
end

sol_table = cell2table(sol_table,"VariableNames",["model_based","T=6","T=10","T=15"]);
sol_table = addvars(sol_table, keywords', 'Before', 1, 'NewVariableNames', 'Design');
h2_table = array2table(h2_table,"VariableNames",["model_based","T=6","T=10","T=15"]);
h2_table = addvars(h2_table, keywords', 'Before', 1, 'NewVariableNames', 'Design');

tv = datetime('now', 'Format','yyyy_MM_dd_HH_mm');
fname = sprintf('H2_epsilon_0.1__%s.mat', tv);
save(".\data\"+fname,  'sol_table', 'h2_table');