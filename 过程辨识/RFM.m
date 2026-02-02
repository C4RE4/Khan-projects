a1 = -1.5;
a2 =  0.7;
b1 =  1.0;
b2 =  0.5;
theta_true = [a1; a2; b1; b2];
N = 1000;          % 数据长度
L = 100;           % 记忆长度
sigma_e = 1;       % 噪声方差
u = randn(N,1);
e = sqrt(sigma_e)*randn(N,1);
y = zeros(N,1);
for k = 3:N
    y(k) = -a1*y(k-1) - a2*y(k-2) ...
           + b1*u(k-1) + b2*u(k-2) + e(k);
end
%% RFM 初始化
theta = 0.001*ones(4,1);
P = 1e7*eye(4);
%% RFM
H = zeros(L,4);
z = zeros(L,1);

for i = 1:L
    k = i + 2;
    H(i,:) = [-y(k-1), -y(k-2), u(k-1), u(k-2)];
    z(i) = y(k);
end

P = inv(H'*H);
theta = P*H'*z;
theta_hist = zeros(4,N);

for k = L+3:N

    %% —— 新数据加入 ——
    h_new = [-y(k-1); -y(k-2); u(k-1); u(k-2)];
    z_new = y(k);

    K = P*h_new / (1 + h_new'*P*h_new);
    theta = theta + K*(z_new - h_new'*theta);
    P = (eye(4) - K*h_new')*P;

    %% —— 去掉最老数据 ——
    h_old = [-y(k-L-1); -y(k-L-2); u(k-L-1); u(k-L-2)];
    denom = 1 - h_old'*P*h_old;
    P = (eye(4) + (P*h_old*h_old')/denom)*P;
    theta = theta - P*h_old*(y(k-L) - h_old'*theta);

    theta_hist(:,k) = theta;
end
%% 输出
theta_est = theta_hist(:,end)

disp('真实参数：')
disp(theta_true)

disp('估计参数：')
disp(theta_est)
