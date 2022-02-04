function E = entropy_base_e(Distr)
% Calculate the Entropy using log based on e, instead 2.
% Input: 
%   Distr: a matrix, with each row is a distribution.
% Output:
%   E: a column vector, corresponding to each the entropy of each distr.

Distr = Distr + (Distr==0);
E = -1 * sum(Distr .* log(Distr), 1); % sum the rows
