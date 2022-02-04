function Data = GenDir(a,n)
% Generat samples from single Dirichlet distribution.
% Input:
%   a: M-by-1 vector. Dirichlet parameter
%   n: Number of samples.
% Output:
%   Data: M-by-N matrix with each column being one sample.
% Note: This function is adapted from Minka's function: dirichelt_sample.m

if size(a,2) ~= 1
    fprintf('$ input a should be a column vector');
    Data=[];
    exit(1);
end

Data = gamrnd(repmat(a, 1, n),1);
s = sum(Data,1);
s(find(s == 0)) = 1;
Data = Data/diag(s); % faster than repmat version.
