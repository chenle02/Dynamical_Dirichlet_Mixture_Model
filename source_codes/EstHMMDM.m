function [A, B, C, Pi] = EstHMMDM(Data, K, M)
% Estimate the parameters of HMM+DM by EM algorithm.
% Input: 
%   Data: 1-by-D cell, with Data{d} is N-by-Td data matrix. 
%         D is the number of sequences;
%         N is data dimension; 
%         Td is number of samples in d^th sequence.
%   K: number of hidden states.
%   M: number of mixture components.
% Output: 
%     A: M-by-N-by-K positive array.
%        A(m,:,k) is the Dirichlet corresponding to (m,k)^th Dirichlet.
%     B: K-by-K probability transition matrix. 
%        B(i,j) = p(h=j|h=i). Then each row of B should sum to 1.
%     C: K-by-M probability matrix.
%        C(i,j) = p(m=j|h=i). Then each row of C should sum to 1. 
%     Pi: K-by-1 column vector, initial probability of hidden states.


% % % --- Debug ------
% M = 3;
% N = 4;
% K = 2;
% D = 6;
% BaseNum = 3000;
% Number = sort(round(rand(D,1)*BaseNum));
% A = round(rand(M,N,K) * 10+1);
% B = rand(K,K);
% for k=1:K; B(k,:) = B(k,:)/sum(B(k,:));end
% C = rand(K,M);
% for k=1:K; C(k,:) = C(k,:)/sum(C(k,:));end
% % Pi = rand(K,1);
% Pi = [0.3;0.7]
% Pi = Pi/sum(Pi);
% Data={};
% Ind = {};
% for d=1:D
%     [Data{d} Ind{d}] = GenDynamicMixtureDir(A,B,C,Pi,Number(d));
% end
% A0=A; B0=B; C0=C; Pi0=Pi;
% save(sprintf('HMMDM_%d',BaseNum));
% % ----------------


tic;fprintf('# Initialize ....');
D = size(Data,2);
N = size(Data{1},1);
AllData = [];
for d=1:D
    AllData = [AllData Data{d}];
end
[A p] = MomentMatchingInitDM(K*M, AllData);
A = permute( reshape(A',[N M K]), [2 1 3]);
Aini = A;
clear AllData;
fprintf('   %6.2f sec. \n',toc);


LogData={};
TotalNumber = 0;
for d=1:D
    LogData{d} = log(Data{d});
    TotalNumber = TotalNumber + size(Data{d},1);
end

Gm={};
Xi={};
A = Aini;
% % A=A0;B=B0;C=C0;Pi=Pi0;
% A = rand(M,N,K)+1;
Pi = rand(K,1);Pi=Pi/sum(Pi);
B = rand(K,K); B = inv(diag(sum(B,2))) * B;
C = rand(K,M); C = inv(diag(sum(B,2))) * C; 

Lnew = 1;
Lold = 0;
nIter=0;
L = [];
Lentr=[];
Lenergy=[];
while abs(Lnew - Lold) > 1e-7
    
    % E step -- Inference for each sequence.
    for d = 1:D
        % -- bottle neck!!! dataLikelihood_DM ---
        obslik = dataLikelihood_DM(A,Data{d});
        [Gm{d}, Xi{d}] = forback(B,C,Pi,obslik);
    end
    
    % M step
    cumPi = zeros(K,1);
    cumXi = zeros(K,K);
    cumGm = zeros(K,M);
    NX = zeros(N,K,M);
    for d=1:D
        for m=1:M
            for k=1:K
                NX(:,k,m) = NX(:,k,m) + LogData{d}*Gm{d}(:,m,k);
            end
        end
        cumPi = cumPi + permute(sum(Gm{d}(1,:,:),2),[3 2 1]);
        cumXi = cumXi + permute(sum(Xi{d},1),[2 3 1]);
        cumGm = cumGm + permute(sum(Gm{d},1),[3,2,1]);
    end
    Pi = cumPi/sum(cumPi);
    B = inv(diag(sum(cumXi,2))) * cumXi;
    C = inv(diag(sum(cumGm,2))) * cumGm;
    for k=1:K
        for m=1:M
            A(m,:,k) = EstDirchlet(A(m,:,k)',NX(:,k,m)/cumGm(k,m))';
        end
    end
    
    % Calculate the low ound.
    Entr = 0;
    for d=1:D
        Td = size(Gm{d},1);
        qq = entropy_base_e(permute(sum(Gm{d},2),[3 1 2]));
        Entr =   Entr   ... 
               + sum(entropy_base_e(reshape(Gm{d},[Td,  K*M,1])'))  ...
               + sum(entropy_base_e(reshape(Xi{d},[Td-1,K*K,1])'))  ...
               - 2 * sum(qq) + qq(1) + qq(Td);
    end
    Energy =   sum(cumPi.*log(Pi))      ...
        + sum(sum(cumGm.*log(C)))  ...
        + sum(sum(cumXi.*log(B)));
    Lenergy = [Lenergy Energy/TotalNumber];
    Energy = Energy + sum(sum(cumGm.*log(gamma(permute(sum(A,2),[3 1 2])))))   ...
        - sum(sum(cumGm.*permute(sum(log(gamma(A)),2),[3 1 2])))   ...
        + sum(sum(sum( (permute(A,[2 3 1])-1) .* NX )));
    Lold = Lnew;
    Lentr = [Lentr Entr/TotalNumber];
    Lnew = (Energy + Entr)/TotalNumber;
    L = [L Lnew];
    if Lold > Lnew
        fprintf('! Warning: decreasing now.\n');
    end
    
    % --- Check maximal iteration ---
    nIter = nIter + 1;
    if nIter >120
        fprintf('! Terminate on maximum %d iterations\n', nIter-1);
        fprintf('# No.  %-d      iter., %6.2f sec.   Average log likelihood	%f \n', [nIter-1, toc, Lnew]);
        break;
    elseif mod(nIter-1,10) == 0
        fprintf('# No.  %-d      iter., %6.2f sec.   Average log likelihood	%f \n', [nIter-1, toc, Lnew]);
    end
end


% % --- Plot Figure (for debug)---------------
% figure;
% set(0,'defaulttextinterpreter','none');
% axes('fontsize',10);
% hold on;
% plot(L,'r-');
% plot(Lentr,'g--'); 
% plot(Lenergy,'m:'); 
% plot(L - Lentr - Lenergy,'k:');
% plot(L - Lentr,'b-.');
% legend('Loglikelihood','Entropy','Energy of Trans.','Energy of DM','Total Energy','location','SouthEast');
% xlabel('Number of EM iteration','fontsize',12);
% ylabel('Energy, Log-likelihood, Entropy','fontsize',12);
% laprint(f,'HMMMD5000_6sequences','width',10,'scalefonts','off');
% --- ----------- ----------------------------

