clc; clear; close all;
%% 参数设置（教材例5.5）
L = 480;                  % 数据长度
theta_true = [-1.5; 0.7; 1.0; 0.5];   % [a1 a2 b1 b2]^T

P = 1e6 * eye(4);         % P(0)
theta_hat = 0.001 * ones(4,1);  % theta(0)

Lambda = 1;               % 加权因子

%% 生成输入信号 u(k)

u = randn(L,1);           % 白噪声输入
v = randn(L,1);           % 测量噪声

z = zeros(L,1);           % 输出初始化

%% 系统仿真（真实系统）

for k = 3:L
    z(k) = ...
        -theta_true(1)*z(k-1) ...
        -theta_true(2)*z(k-2) ...
        +theta_true(3)*u(k-1) ...
        +theta_true(4)*u(k-2) ...
        +v(k);
end

%% RLS 递推最小二乘

theta_history = zeros(4,L);

for k = 3:L
    % 回归向量 h(k)
    h = [-z(k-1);
         -z(k-2);
          u(k-1);
          u(k-2)];
    
    % 增益矩阵 K(k)
    K = P*h / (h'*P*h + 1/Lambda);
    
    % 预测误差
    e = z(k) - h'*theta_hat;
    
    % 参数更新
    theta_hat = theta_hat + K*e;
    
    % 协方差更新
    P = (eye(4) - K*h')*P;
    
    % 保存历史
    theta_history(:,k) = theta_hat;
end

%% 结果绘图

figure;
subplot(2,2,1)
plot(theta_history(1,:), 'LineWidth',1.5); hold on;
yline(theta_true(1),'r--'); grid on;
title('a_1 估计')

subplot(2,2,2)
plot(theta_history(2,:), 'LineWidth',1.5); hold on;
yline(theta_true(2),'r--'); grid on;
title('a_2 估计')

subplot(2,2,3)
plot(theta_history(3,:), 'LineWidth',1.5); hold on;
yline(theta_true(3),'r--'); grid on;
title('b_1 估计')

subplot(2,2,4)
plot(theta_history(4,:), 'LineWidth',1.5); hold on;
yline(theta_true(4),'r--'); grid on;
title('b_2 估计')

sgtitle('RLS 参数辨识（教材例 5.5）')

%% ======================
% 6. 最终结果输出
%% ======================
disp('真实参数：')
disp(theta_true.')

disp('RLS 估计参数：')
disp(theta_hat.')
