function v = compute_residual(theta,z,u,n)

L = length(z);

a = theta(1:n);
b = theta(n+1:2*n);
d = theta(2*n+1:3*n);

v = zeros(L,1);

for k = n+1:L
    
    Az = 0;
    Bu = 0;
    Dv = 0;
    
    for i = 1:n
        Az = Az + a(i)*z(k-i);
        Bu = Bu + b(i)*u(k-i);
        Dv = Dv + d(i)*v(k-i);
    end
    
    v(k) = z(k) + Az - Bu - Dv;
end

end
