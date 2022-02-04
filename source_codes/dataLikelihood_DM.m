function obslik = dataLikelihood_DM(A,data,isLog)
% Calculate the data likelihood for the Dirichlet mixture.
% Input:
%     A: M-by-N-by-K matrix, parameters of DM. 
%        M*K is number of mixture components.
%        N is sample dimension.
%        When K =1, A is a matrix; When K=1 and M=1, it is a single Dir.
%     data: N-by-T matrix. T is the sample number.
%     isLog: 0 - output likelihood (Default), otherwise log likelihood.
% Output:
%     obslik: T-by-M-by-K matrix.
%             obslik(t,m,k) is the t^th sample's (log-)likelihood on
%             (m,k)^th mixture components.


% %--- Debug -----
% M = 5;
% N = 8;
% K = 7;
% T = 12;
% A = rand(M,N,K) * 10;
% data = rand(N,T);
% for t=1:T; data(:,t) = data(:,t)/sum(data(:,t));end
% %---------------

if nargin < 3
    isLog = 0;
end

[M,N,K] = size(A);
T = size(data,2);
obslik = zeros(T,M,K);

for k=1:K
    for m=1:M
        obslik(:,m,k) = Dirichlet_loglike(A(m,:,k)',data)';
        if isLog == 0
            obslik(:,m,k) = exp(obslik(:,m,k));
        end
    end
end