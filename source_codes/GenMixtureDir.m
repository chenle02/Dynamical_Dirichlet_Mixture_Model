function Data = GenMixtureDir(A,Pi,Number)
% Generate samples from Dirichlet mixture model.
% Input:
%   A:  M-by-N matrix. Parameters for Mixture of Dirichlet
%       each row is one Dirichlet.
%       N: Data dimension
%       M: number of mixture
%   Number: Number of samples to generate.
%   Pi: M-by-1 vector. Prior distribution for each Dirichlet.
% Output:
%   Data: N-by-Number matrix with each column being one sample. 


% --- Debug -----
% M = 3;
% N = 5;
% Pi = rand(M,1);Pi=Pi/sum(Pi);
% A = round(rand(M,N)*10) + 1;
% Number = 1000;
% ---------------

M = size(A,1);
K = round(Pi*Number);
Data=[];
for m=1:M
    Data = [Data GenDir(A(m,:)',K(m))];
end

