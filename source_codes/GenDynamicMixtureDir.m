function [Data Ind]= GenDynamicMixtureDir(A,B,C,Pi,Number)
% Generate one sequence of random samples from HMM+DM model.
% Input:
%     A: M-by-N-by-K positive array.
%        M is number of Dirichlet components.
%        N is sample's dimension.
%        K is the number of hidden states.
%        A(m,:,k) is the Dirichlet corresponding to (m,k)^th Dirichlet.
%     B: K-by-K probability transition matrix. 
%        B(i,j) = p(h=j|h=i). Then each row of B should sum to 1.
%     C: K-by-M probability matrix.
%        C(i,j) = p(m=j|h=i). Then each row of C should sum to 1. 
%     Pi: K-by-1 column vector, initial probability of hidden states.
%     Number: Number of samples to generate.
%             T is the sample number.
%             obslik(t,m,k) is t^th samples likelihood on (m,k)^th component. 
% Output:
%   Data: N-by-Number matrix with each column is one sample. 
%   Ind:  2-by-Number matrix, with Ind(1,t) in {1,...,K} denoting which
%         hidden state sample t belongs. And Ind(2,t) in {1,...,M}
%         denoting which Dirichlet generating sample t.


% --- Debug ---------
% M = 3;
% N = 3;
% K = 2;
% Number = 2000;
% A = round(rand(M,N,K) * 10+1);
% B = rand(K,K);
% for k=1:K; B(k,:) = B(k,:)/sum(B(k,:));end
% C = rand(K,M);
% for k=1:K; C(k,:) = C(k,:)/sum(C(k,:));end
% % Pi = rand(K,1);
% Pi = [0.3;0.7];
% Pi = Pi/sum(Pi);
% -------------------

[M,N,K] = size(A);

cumPi = [0; cumsum(Pi)];
cumB = [repmat(0,K,1) cumsum(B,2)];
cumC = [repmat(0,K,1) cumsum(C,2)];
Data = zeros(N,Number);
Ind  = zeros(2,Number);
for t=1:Number
    r = rand(1);    
    if t == 1
        Ik = find(cumPi >= r)-1;
    else
        Ik = find(cumB(Ik(1),:) >= r)-1;
    end
    r = rand(1);    Im = find(cumC(Ik(1),:) >= r)-1;
    Data(:,t) = GenDir(A(Im(1),:,Ik(1))',1);
    Ind(:,t) = [Ik(1);Im(1)];
end
