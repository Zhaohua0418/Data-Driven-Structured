LastName = {'Sanchez';'Johnson';'Li';'Diaz';'Brown'};
Age = [38;43;38;40;49];
Smoker = logical([1;0;1;0;1]);
Height = [71;69;64;67;64];
Weight = [176;163;131;133;119];
BloodPressure = [124 93; 109 77; 125 83; 117 75; 122 80];

T = table(LastName,Age,Smoker,Height,Weight,BloodPressure)

% 示例表格
T = table([1; 2; 3], {'A'; 'B'; 'C'}, 'VariableNames', {'Number', 'Letter'});
% 新列数据（确保与表格行数一致）
newColumn = [10; 20; 30];

% 将新列添加到表格最左侧
T = addvars(T, newColumn, 'Before', 1, 'NewVariableNames', 'NewColumn');

