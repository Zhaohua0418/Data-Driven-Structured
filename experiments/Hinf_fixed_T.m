clear all
rng(45, 'twister');
nx = 3;
nu = 2;
ny = 3;
nd = 2;

S = [1 1 0;...
     1 1 0];


% epsilon list
epsilon_list = [0.01; 0.05; 0.1; 0.15; 0.2; 0.5];


T = 50;


% system parameters
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

SIM_list = cell(1,length(epsilon_list)); 
SIM = simulator(sys);
for j=1:length(epsilon_list)
    SIM_list{j} = SIM.run(epsilon_list(j),T);
end

% keywords
keywords = ["diag_unstructured","my_unstructured", "diag_structured", "my_structured_energy_bound","diag_structured_instantaneous_bound","my_structured"];
classes = {Hinf_module_diag_unstructured(sys),Hinf_module_my_unstructured(sys),Hinf_module_diag_structured(sys,S),Hinf_module_my_structured_energy_bound(sys,S),Hinf_module_diag_structured_instantaneous_bound(sys,S),Hinf_module_my_structured(sys,S)};

d = dictionary(keywords,classes);
sol_table = cell(length(keywords),1+length(epsilon_list)); 
hinf_table = NaN*zeros(length(keywords),1+length(epsilon_list));

for i=1:length(keywords)
    sol_table{i,1} = classes{i}.model_based_Hinf(); % cell一定要用{}
    hinf_table(i,1) = sol_table{i,1}.hinf;
    for j=1:length(epsilon_list)
        SIM_out = SIM_list{j};
        try
            sol_table{i,j+1} = classes{i}.data_driven_Hinf(SIM_out);
            hinf_table(i,j+1) = sol_table{i,j+1}.hinf;
        catch
            sol_table{i,j+1} = NaN;
        end
    end
end

sol_table = cell2table(sol_table,"VariableNames",["model_based","eps=0.01","eps=0.05","eps=0.1","eps=0.15","eps=0.2","eps=0.5"]);
sol_table = addvars(sol_table, keywords', 'Before', 1, 'NewVariableNames', 'Design');
hinf_table = array2table(hinf_table,"VariableNames",["model_based","eps=0.01","eps=0.05","eps=0.1","eps=0.15","eps=0.2","eps=0.5"]);
hinf_table = addvars(hinf_table, keywords', 'Before', 1, 'NewVariableNames', 'Design');
% 
tv = datetime('now', 'Format','yyyy_MM_dd_HH_mm');
fname = sprintf('Hinf_T_%d__%s.mat', T,tv);
save(".\data\"+fname,  'sol_table', 'hinf_table');

