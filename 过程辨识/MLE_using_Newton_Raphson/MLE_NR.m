function theta = MLE_NR(y, u, MaxIter)
% -------------------------------------------------
% Maximum Likelihood Estimation with Newton-Raphson
% ARX structure (Example 7.5)
% -------------------------------------------------

N = length(y);
theta = zeros(4,1);     % 初值（非常关键：必须非 NaN）

for iter = 1:MaxIter

    grad = zeros(4,1);
    H = zeros(4,4);

    for k = 3:N

        % 回归向量
        phi = [-y(k-1);
               -y(k-2);
                u(k-1);
                u(k-2)];

        % 预测误差
        eps = y(k) - phi.'*theta;

        % 梯度
        grad = grad - phi*eps;

        % Hessian 近似
        H = H + phi*phi.';
    end

    % 正则化
    H = H + 1e-6*eye(4);

    % Newton 更新
    delta = H \ grad;
    theta = theta - delta;
end
end

