function [A, Pi] = MomentMatchingInitDM(M,DataSet)
% Initialize the mixture of Dirichlet by Kmeans + Moment Matching.
% Input:
%   M:  Number of mixture components (M >=1).
%   DataSet: N-by-T sample matrix with 
%            N is the data dimension
%            T is the number of samples.
% Output: 
%   A:  M-by-N matrix, with each row corresponding one mixture component.
%   Pi: M-by-1 vector, the prior probabilities for each Dirichlet.

% % --- Debug ----
% M = 3;
% N = 3;
% T = 300;
% DataSet = rand(N,T);
% for t=1:T; DataSet(:,t) = DataSet(:,t)/sum(DataSet(:,t));end

% % --------------


[N,T] = size(DataSet);
% [mu, Sigma, prior] = mixgauss_em(DataSet, M);
% A = zeros(M,N);
% for m=1:M
%     P = inv(Sigma(:,:,m));
%     detS = det(2*pi*Sigma(:,:,m));
%     for t=1:T
%         p(m,t) = exp(-DataSet(:,t)'*P*DataSet(:,t)/2)/detS;
%     end
% end
if M > 1
    Ind = kmeans(DataSet',M);
    p = zeros(M,T);
    for t=1:T
        p(Ind(t),t) = 1;
    end
else
    p = repmat(1,1,T);
end
Pi = sum(p,2);
Pi = Pi/sum(Pi);

%---- Main Loop ----
A = zeros(M,N);
for m=1:M
    y = zeros(1,N);
    z = zeros(1,N);
    s = 0;
    for t=1:T
        s = s + p(m,t);
        y = y + DataSet(:,t)' * p(m,t);
        z = z + DataSet(:,t)'.* DataSet(:,t)' * p(m,t);
    end
    y = y/s;
    z = z/s;
    w = sum(log((y-z)./(z-y.*y)))/N;
    A(m,:) = y * exp(w);
end
