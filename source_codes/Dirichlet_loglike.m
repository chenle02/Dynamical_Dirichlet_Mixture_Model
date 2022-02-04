function p = Dirichlet_loglike(a, data)
% Evaluate the data's log-likelihoods relative to a Dirichlet model.
% Input: 
%   a:    a N-by-1 column vector, Dirichlet parameter.
%   data: a N-by-T matrix, with each column being one sample (sum to one).
% Output:
%   p: a 1-by-T row vector for log likelihoods.
% Note: This procedure is adapted from Minka's dirichlet_logProb.m


% Check parameters
if size(a,2) ~= 1 || size(a,1) ~= size(data,1)
    fprintf('$ Wrong input...Please try again.');
    exit(1);
end

w = warning;
warning('off')
p = (a-1)' * log(data);
p = p' + gammaln(sum(a)) - sum(gammaln(a));
warning(w);
