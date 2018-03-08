function [Phi, Omega, b] = reg_dmd(x,TimeSpan,r)
    XP = x(:,2:end);
    X = x(:,1:end-1);

    [U,S,V] = svd(X,'econ');
    
    r = min(r,size(U,2));
    
    U_r = U(:,1:r);
    V_r = V(:,1:r);
    S_r = S(1:r,1:r);

    A_tilde = U_r' * XP * V_r * (S_r^(-1));

    [W_r,D] = eig(A_tilde);
    lambda = diag(D);

    Phi = XP * V_r * (S_r^(-1)) * W_r;
    
    Omega = log(lambda)/(TimeSpan(2) - TimeSpan(1));

    b = pinv(Phi)*x(:,1); 
end