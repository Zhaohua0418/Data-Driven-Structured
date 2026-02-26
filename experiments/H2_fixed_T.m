clear all
rng(45, 'twister');
cvx_solver mosek
nx = 3;
nu = 2;
ny = 5;
nd = 3;

S = [1 1 0; 
    0 1 1];



% epsilon list
epsilon_list = [0.05; 0.1; 0.15; 0.2];

% T = 6;
% T = 10;
T = 20;
load("H2_info_struct_epsilon_list_T_20.mat");

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

SIM_list = cell(1,length(epsilon_list)); % 目前还没有用到这个，后面这个要自己设计

% keywords
keywords = ["jared_unstructured_energy_bound","my_unstructured_instantaneous_bound", "jared_structured_energy_bound", "my_structured_energy_bound", "jared_structured_instantaneous_bound","my_structured_instantaneous_bound"];
classes = {H2_module_jared_unstructured(sys),H2_module_my_unstructured(sys),H2_module_jared_structured(sys, S),H2_module_my_structured_energy_bound(sys, S),H2_module_jared_structured_instantaneous_bound(sys, S),H2_module_my_structured(sys, S)};
d = dictionary(keywords,classes);
sol_table = cell(length(keywords),1+length(epsilon_list)); 
h2_table = NaN*zeros(length(keywords),1+length(epsilon_list));
for i=1:length(keywords)
    sol_table{i,1} = classes{i}.model_based_H2(); % cell一定要用{}
    h2_table(i,1) = sol_table{i,1}.h2;
    for j=1:length(epsilon_list)
        try
            sol_table{i,j+1} = classes{i}.data_driven_H2(traj{j});
            h2_table(i,j+1) = sol_table{i,j+1}.h2;
        catch % catch infeasibility
            sol_table{i,j+1} = NaN;
        end
    end
end

sol_table = cell2table(sol_table,"VariableNames",["model_based","eps=0.05","eps=0.1","eps=0.15","eps=0.2"]);
sol_table = addvars(sol_table, keywords', 'Before', 1, 'NewVariableNames', 'Design');
h2_table = array2table(h2_table,"VariableNames",["model_based","eps=0.05","eps=0.1","eps=0.15","eps=0.2"]);
h2_table = addvars(h2_table, keywords', 'Before', 1, 'NewVariableNames', 'Design');

tv = datetime('now', 'Format','yyyy_MM_dd_HH_mm');
fname = sprintf('H2_T_%d__%s.mat', T,tv);
save(".\data\"+fname,  'sol_table', 'h2_table');