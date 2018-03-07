function [Phi, Omega, b] = reg_dmd(x,TimeSpan)
   XP = x(:,2:end);
    X = x(:,1:end-1);

    [U,S,V] = svd(X,'econ');

    A_tilde = U' * XP * V * (S^(-1));

    [W,Lambda] = eig(A_tilde);

    Phi = XP * V * (S^(-1)) * W;
    
    Omega = log(diag(Lambda))/(TimeSpan(2) - TimeSpan(1));

    b = pinv(Phi)*x(:,1); 
end