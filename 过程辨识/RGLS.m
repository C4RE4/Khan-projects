%% 参数与初始化
na = 2; nb = 2; nc = 2;
N  = 480;

theta_hat  = zeros(na+nb,1);
theta_c_hat = zeros(nc,1);

Pf = 1e5*eye(na+nb);
Pc = eye(nc);

e_hat = zeros(N,1);

%% 主循环
for k = max([na,nb,nc])+1:N

    % 回归向量
    h = [-z(k-1:-1:k-na);
          u(k-1:-1:k-nb)];

    % 系统参数 RLS  
    Kf = Pf*h/(1 + h'*Pf*h);
    theta_hat = theta_hat + Kf*(z(k) - h'*theta_hat);
    Pf = (eye(length(theta_hat)) - Kf*h')*Pf;

    % 残差
    e_hat(k) = z(k) - h'*theta_hat;

    % 噪声模型 RLS 
    hc = -e_hat(k-1:-1:k-nc);
    Kc = Pc*hc/(1 + hc'*Pc*hc);
    theta_c_hat = theta_c_hat + Kc*(e_hat(k) - hc'*theta_c_hat);
    Pc = (eye(nc) - Kc*hc')*Pc;
end

%% 画图
plot(theta_ls_hist'); grid on;
legend('a1','a2','b1','b2')


