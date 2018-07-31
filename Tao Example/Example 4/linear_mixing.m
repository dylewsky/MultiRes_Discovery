clear variables; close all; clc

rng(111); %seed
load raw_data_3_hiRes.mat;
n = size(x,2);

% generate a random unitary matrix M
X = rand(n)/sqrt(2);
[Q,R] = qr(X);
R = diag(diag(R)./abs(diag(R)));
M = real(Q*R);

x = x * M;
save('raw_data_4_hiRes.mat','x','TimeSpan','M');