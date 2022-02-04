function [A, Pi] = EstMixDirichlet(Data, M)
% Estimate the parameters of the mixture of Dirichlet by EM algorithm.
% Input: 
% 	Data: N-by-T data matrix. N is data dimension; T is number of samples.
%   M: number of mixture components.
% Output: 
%   A:  M-by-N matrix, with each row corresponding to one Dirichlet.
%   Pi: M-by-1 vector, the prior probabilities for each Dirichlet.


% % --- Debug -----
% M = 3;
% N = 5;
% Pi0 = [0.2 0.3 0.5];
% Pi0 = rand(M,1);Pi0=Pi0/sum(Pi0);
% A0 = round(rand(M,N)*10) + 1;
% Number = 2000;
% Data = GenMixtureDir(A0,Pi0,Number);
% % ---------------

%---------- Exp in the paper -----
% M = 3;
% N = 5;
% A0 = [3 3 4 6 5
%        10 7 1 9 10
%        2 6 2 9 10];
% Pi0 = [0.20
%         0.30
%         0.50];
% Number = 2000;
% Data = GenMixtureDir(A0,Pi0,Number);
% save(sprintf('DM_%d',Number));
% --------------------------------


[A, Pi] = MomentMatchingInitDM(M,Data);
Pi = rand(M,1);Pi=Pi/sum(Pi);
LogData = log(Data);
Lnew = 1;
Lold = 0;
L = [];
Lentr=[]
while abs(Lnew - Lold) > 1e-8
    % E step
    obslik = dataLikelihood_DM(A,Data,0) * diag(Pi);
    for t = 1:size(obslik,1)
        obslik(t,:) = obslik(t,:) / sum(obslik(t,:));
    end
    
    % M step
    g = sum(obslik)'; Pi = g/sum(g);
    nx = LogData * obslik / diag(g);
    for m=1:M
        A(m,:) = EstDirchlet(A(m,:)',nx(:,m))';
    end

    % Calculate the low ound.
    Entr = sum(entropy_base_e(obslik'));
    Energy = g'*( log(Pi) +  log(gamma(sum(A,2))) - sum(log(gamma(A)),2) + sum((A-1).*nx',2) );
    Lold = Lnew;
    Lnew = Energy + Entr;
    L = [L Lnew];
end

% ---- Graphical show ------
f = figure;
set(0,'defaulttextinterpreter','none')
axes('fontsize',10);
plot(L/Number); hold on;
ylabel('Aver. Data log likelihood','fontsize',12);
xlabel('Number of EM Iterations','fontsize',12);
% Title(sprintf('%d %d-Dim. Samples, %d mixture components',[Number N M]),'fontsize',12);
% laprint(f,'StaticMD2000Samples','width',7,'scalefonts','off');
% ---------------------------

