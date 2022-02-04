% Script to test EstHMMDM.m

% -- Generate data from HMM+DM model first
fprintf('# Generating data first ....\n');
M = 3;  % number of Dirichlets.
N = 10; % data dimension
K = 5;  % number of hidden states h
D = 6;  % number of sequences
BaseNum = 60;
Number = sort(round(rand(D,1)*BaseNum));
A = round(rand(M,N,K) * 10+1);
B = rand(K,K);
for k=1:K; B(k,:) = B(k,:)/sum(B(k,:));end
C = rand(K,M);
for k=1:K; C(k,:) = C(k,:)/sum(C(k,:));end
Pi = rand(K,1);
% Pi = [0.3;0.7]
Pi = Pi/sum(Pi);
Data={};
Ind = {};
for d=1:D
    [Data{d} Ind{d}] = GenDynamicMixtureDir(A,B,C,Pi,Number(d));
end
A0=A; B0=B; C0=C; Pi0=Pi;
save(sprintf('HMMDM_%d',BaseNum));


% -- Estimating....
fprintf('# Estimating now....\n');
[A, B, C, Pi] = EstHMMDM(Data, K, M);
