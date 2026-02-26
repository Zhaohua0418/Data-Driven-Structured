clear all
rng(45, 'twister');
cvx_solver mosek
nx = 3;
nu = 2;
ny = 3;
nd = 2;

S = [1 1 0;...
     1 1 0];

% epsilon list
% epsilon_list = [0.05];

T = 50;
% t = 10;

T_list = [10 20 40];

fname = sprintf('Hinf_info_struct_epsilon_0.05_T_%d.mat', T);
load(fname);


Hinf_info_struct_copy = Hinf_info_struct;


sys = struct(...
    'A', [0.8 0.2  0.1;...
          0.1 0.7 -0.3;...
         -0.3 0.5  0.9],...
    'B', [1 0;...
          0 1;...
          1 1],...
    'G', [0.3 0.1;...
          0.2 0.2;...
          0.1 0.3],...
    'C', [1 0 0;...
          0 1 0;...
          0 0 1],...
    'D', [0.1 0.2;...
          0.3 0.1;...
          0.2 0.1],...
    'H', [0.1 0.1;...
          0.2 0.2;...
          0.3 0.3],...
    'nx',nx,'nu',nu,'ny',ny,'nd',nd);





%%
hinf_module = Hinf_module_my_structured(sys,S);

%% model-based control
MBHinf_sol = hinf_module.model_based_Hinf();
MBhinf = MBHinf_sol.hinf;

%% data-driven control
keywords = ["diag_unstructured","my_unstructured", "diag_structured", "my_structured_energy_bound","diag_structured_instantaneous_bound","my_structured"];
classes = {Hinf_module_diag_unstructured(sys),Hinf_module_my_unstructured(sys),Hinf_module_diag_structured(sys,S),Hinf_module_my_structured_energy_bound(sys,S),Hinf_module_diag_structured_instantaneous_bound(sys,S),Hinf_module_my_structured(sys,S)};
d = dictionary(keywords,classes);
sol_table = cell(length(keywords),1+length(T_list)); 
hinf_table = NaN*zeros(length(keywords),1+length(T_list));

for i=1:length(keywords)
    sol_table{i,1} = classes{i}.model_based_Hinf(); % cell一定要用{}
    hinf_table(i,1) = sol_table{i,1}.hinf;
    for j=1:length(T_list)
        Hinf_info_struct_copy.Xp = Hinf_info_struct.Xp(:,1:T_list(j));
        Hinf_info_struct_copy.Xm = Hinf_info_struct.Xm(:,1:T_list(j));
        Hinf_info_struct_copy.Um = Hinf_info_struct.Um(:,1:T_list(j));
        Hinf_info_struct_copy.T = T_list(j);
        try
            sol_table{i,j+1} = classes{i}.data_driven_Hinf(Hinf_info_struct_copy);
            hinf_table(i,j+1) = sol_table{i,j+1}.hinf;
        catch % catch infeasibility
            sol_table{i,j+1} = NaN;
        end
    end
end

sol_table = cell2table(sol_table,"VariableNames",["model_based","T=10","T=20","T=40"]);
sol_table = addvars(sol_table, keywords', 'Before', 1, 'NewVariableNames', 'Design');
hinf_table = array2table(hinf_table,"VariableNames",["model_based","T=10","T=20","T=40"]);
hinf_table = addvars(hinf_table, keywords', 'Before', 1, 'NewVariableNames', 'Design');

tv = datetime('now', 'Format','yyyy_MM_dd_HH_mm');
fname = sprintf('Hinf_epsilon_0.05__%s.mat', tv);
save(".\data\"+fname,  'sol_table', 'hinf_table');