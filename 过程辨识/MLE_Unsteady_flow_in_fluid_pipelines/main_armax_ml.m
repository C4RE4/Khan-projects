%% 真实系统参数（可修改）
n = 2;                      % 阶次
a_true = [ -1.5  0.7 ];     % A 多项式参数
b_true = [  1.0  0.5 ];     % B 多项式参数
d_true = [  0.3 -0.2 ];     % D 多项式参数

L = 1000;                   % 数据长度
sigma_v = 0.2;              % 白噪声标准差

%% 输入信号
u = randn(L,1);             % 高斯白噪声输入

%% ARMAX系统输出
v = sigma_v * randn(L,1);   % 白噪声

z = zeros(L,1);

for k = n+1:L
    
    Az = 0;
    Bu = 0;
    Dv = 0;
    
    for i = 1:n
        Az = Az + a_true(i)*z(k-i);
        Bu = Bu + b_true(i)*u(k-i);
        Dv = Dv + d_true(i)*v(k-i);
    end
    
    z(k) = -Az + Bu + v(k) + Dv;
end

%% 极大似然估计
% 初始参数
theta0 = zeros(3*n,1);

options = optimoptions('fminunc',...
    'Algorithm','quasi-newton',...
    'Display','iter',...
    'MaxIterations',200);

theta_hat = fminunc(@(theta) cost_function(theta,z,u,n), theta0, options);

%% 输出


a_hat = theta_hat(1:n)';
b_hat = theta_hat(n+1:2*n)';
d_hat = theta_hat(2*n+1:3*n)';

disp('真实参数：')
disp([a_true; b_true; d_true])

disp('估计参数：')
disp([a_hat; b_hat; d_hat])

%% 残差分析

v_est = compute_residual(theta_hat,z,u,n);

figure;
plot(v_est)
title('Estimated Residual')

figure;
autocorr(v_est)
title('Residual Autocorrelation')
