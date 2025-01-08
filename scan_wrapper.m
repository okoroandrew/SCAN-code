function [Omega, pred] = scan_wrapper(data) 

[t1,f] = size(data);
t = 2^floor(log2(t1));
data = data(1:t, :);

addpath(genpath("func"))
addpath(genpath("precomputedH"))

%% use precomputed H
if f == 512
    pre_computed_H = load("H_512b.mat");
    H_all = pre_computed_H.H;
    HatomInfo = pre_computed_H.atomInfo;

elseif f == 256
    pre_computed_H = load("H_256.mat");
    H_all = pre_computed_H.H;
    HatomInfo = pre_computed_H.atomInfo;

elseif f == 1024
    pre_computed_H = load("H_1024.mat");
    H_all = pre_computed_H.H;
    HatomInfo = pre_computed_H.atomInfo;

elseif f == 768
    pre_computed_H = load("H_768b.mat");
    H_all = pre_computed_H.H;
    HatomInfo = pre_computed_H.atomInfo;

else
    [H, atomInfo] = H_sparse_gen(f);
    H_var = strcat("H_", num2str(f), "b.mat");
    save(H_var, "H", "atomInfo", '-v7.3');
    
    H_all = H;
    HatomInfo = atomInfo;
end


%% encoding params
lambda = 0.0000001; 
rho = 10*lambda;

threshold=0.15;

%% generate a dictionary
Phi = generate_haar(t);
Phi = Phi';

beta = 1;


% Call scan main
[Omega,pred] = scan_samp(data,H_all,HatomInfo,Phi,threshold,rho,lambda, beta);







