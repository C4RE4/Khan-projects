%  Example 7.5  数据生成
N = 300;                 
u = idinput(N,'prbs');   

% 有色噪声 
e = randn(N,1);
v = zeros(N,1);
for k = 2:N
    v(k) = 0.6*v(k-1) + e(k);
end

% 真参数（仅用于对比）
theta_true = [-1.5; 0.7; 1.0; 0.5];

% 系统输出
y = zeros(N,1);
for k = 3:N
    y(k) = ...
        -theta_true(1)*y(k-1) ...
        -theta_true(2)*y(k-2) ...
        +theta_true(3)*u(k-1) ...
        +theta_true(4)*u(k-2) ...
        + v(k);
end


%  MLE + Newton-Raphson

MaxIter = 50;                 % 最大迭代次数
theta_hat = MLE_NR(y, u, MaxIter);


%  结果输出


disp('True parameters:')
disp(theta_true.')

disp('Estimated parameters:')
disp(theta_hat.')
