clear variables; close all; clc

rng(111); %seed
load raw_data_5_unmixed.mat;

uv_ratio = 1; %ratio of u and v in linear combination

n = size(uv,2);
m = size(uv,1);

nVars_out = 4; %dimension of space to map into

A = randn(n,nVars_out);

x = uv * A;

save('raw_data_5.mat','x','TimeSpan','A');