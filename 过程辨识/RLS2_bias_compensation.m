%% 系统参数
N = 480;
u = idinput(N,'prbs',[0 1],[-1 1]);

% 真参数
a1 = -1.5; a2 = 0.7;
b1 = 1.0;  b2 = 0.5;

% 有色噪声
v = randn(N,1);
Nf = [1 0.5];
Df = [1 -0.8];
e = filter(Nf,Df,v);

% 输出
z = zeros(N,1);
for k = 3:N
    z(k) = -a1*z(k-1)-a2*z(k-2) ...
           + b1*u(k-1)+b2*u(k-2) ...
           + e(k);
end

%% RLS + 偏差补偿
na = 2; nb = 2;
theta_ls = zeros(4,1);
theta_c  = zeros(4,1);

P = 1e5*eye(4);
J = 0;

D = diag([1 1 0 0]);

theta_ls_hist = zeros(4,N);
theta_c_hist  = zeros(4,N);

for k = 3:N
    h = [-z(k-1); -z(k-2); u(k-1); u(k-2)];
    
    % RLS
    K = P*h/(1+h'*P*h);
    err = z(k) - h'*theta_ls;
    theta_ls = theta_ls + K*err;
    P = (eye(4)-K*h')*P;
    
    % J(k)
    J = J + err^2/(1+h'*P*h);
    
    % 噪声方差估计
    sigma2 = J/(k*(1+theta_c'*D*theta_ls));
    
    % 偏差补偿
    theta_c = theta_ls + k*sigma2*P*D*theta_c;
    
    % 参数随时间的收敛轨迹
    theta_ls_hist(:,k) = theta_ls;
    theta_c_hist(:,k)  = theta_c;
end

%% 绘图
figure;
subplot(2,1,1)
plot(theta_ls_hist')
title('LS 参数估计')

subplot(2,1,2)
plot(theta_c_hist')
title('偏差补偿 LS 参数估计')

