% 定义数据
t = 10:5:200;

% 创建图形
figure;
plot(t, my_running_time(2:40), 'b-', 'LineWidth', 2); % 蓝色实线，加粗
hold on;
plot(t, jared_running_time(2:40), 'r--', 'LineWidth', 2); % 红色虚线，加粗

% 设置横轴和纵轴标签
xlabel('$\textit{T}$', 'Interpreter', 'latex');
ylabel('Computation time (s)');

% 添加图例
legend('Our Running Time', 'Running Time in [15]', 'Location', 'best');

% 设置图形其他属性（可选，增强美观性）
grid on; % 添加网格线
set(gca, 'FontSize', 12); % 设置字体大小
title('Computation Time Comparison'); % 添加标题（可选）

hold off;