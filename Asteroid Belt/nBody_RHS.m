function dzdt = nBody_RHS(z,m,G)
    n = length(z)/4;
    x = [z(1:n) z(n+1:2*n)];
    v = [z(2*n+1:3*n) z(3*n+1:4*n)];
    f = zeros(n,2);
    
    for j = 1:n
        for k = 1:j
            if j == k
                continue
            end
            delta_x = x(j,:) - x(k,:);
            f_jk = G*m(j)*m(k)*delta_x / norm(delta_x)^3;
            f(j,:) = f(j,:) - f_jk;
            f(k,:) = f(k,:) + f_jk;
        end
    end
    
    dxdt = v;
    dvdt = f./repmat(m,1,2);
    dzdt = [dxdt(:,1); dxdt(:,2); dvdt(:,1); dvdt(:,2)];
end