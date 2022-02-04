function [Gm,Xi] = forback(B,C,Pi,obslik)
% Inference procedure for general HMM + Mixture of density model.
% Input:
%     B: K-by-K probability transition matrix. 
%        K is the number of hidden states.
%        B(i,j) = p(h=j|h=i). Then each row of B should sum to 1.
%     C: K-by-M probability matrix.
%        M is the number of mixture densities. 
%        C(i,j) = p(m=j|h=i). Then each row of C should sum to 1. 
%     Pi: K-by-1 column vector, initial probability of hidden states.
%     obslik: T-by-M-by-K likelihood (not log likelihood) arrays.
%             T is the sample number.
%             obslik(t,m,k) is t^th samples likelihood on (m,k)^th component. 
% Output: the smoothed states
%     Gm: T-by-M-by-K matrix. Gm(t,m,k) = p(h_t = k,m_t = m |X_1,...,X_T),
%     Xi: (T-1)-by-K-by-K matrix. Xi(t,k1,k2)=p(h_t=k1,h_{t+1}=k2 |X_1,...,X_T), 
%          t=1,...,T-1.


% % -- Debug ------
% M = 3;
% N = 3;
% K = 2;
% T = 10;
% A = round(rand(M,N,K) * 10 +1);
% % % % data = rand(N,T);
% % % % for t=1:T; data(:,t) = data(:,t)/sum(data(:,t));end
% B = rand(K,K);
% for k=1:K; B(k,:) = B(k,:)/sum(B(k,:));end
% C = rand(K,M);
% for k=1:K; C(k,:) = C(k,:)/sum(C(k,:));end
% Pi = rand(K,1);
% Pi = Pi/sum(Pi);
% [data Ind]= GenDynamicMixtureDir(A,B,C,Pi,T)
% obslik = dataLikelihood_DM(A,data);
% % ---------------

[K,M] = size(C);
T = size(obslik,1);

% ---- Forward pass: Normalized Alpha -------
Alpha = zeros(T,M,K);
for t=1:T
    for m=1:M
        for k =1:K
            if t == 1
                Alpha(t,m,k) = Pi(k) * C(k,m) * obslik(t,m,k);
            else
                for m1 = 1:M
                    for k1 =1:K
                        Alpha(t,m,k) = Alpha(t,m,k) + Alpha(t-1,m1,k1) * B(k1,k) * C(k,m) * obslik(t,m,k);
                    end
                end
            end
        end
    end
    % Normalize 
    Alpha(t,:,:) = Alpha(t,:,:) / sum(sum(Alpha(t,:,:),3),2);
end

% ---- Backward pass: Xi, Gm -------
Gm = Alpha;
Xi = zeros(T-1,K,K);
for t = T-1:-1:1
    % --- Normalization constatnt ----
    for k1=1:K
        NC(k1) = 0;
        for k =1:K
            for m=1:M
                NC(k1) = NC(k1) + Alpha(t,m,k) * B(k,k1);
            end
        end
    end
    
    % ----- Main loop ---------
    for k=1:K
        for m=1:M
            Gm(t,m,k) = 0;
            for k1 = 1:K
                for m1=1:M
                    Gm(t,m,k) = Gm(t,m,k) + Alpha(t,m,k) * B(k,k1) * Gm(t+1,m1,k1) / NC(k1);
                end
            end
        end
    end
    for k=1:K
        for k1=1:K
            Xi(t,k,k1) = 0;
            for m=1:M
                for m1=1:M
                    Xi(t,k,k1) = Xi(t,k,k1) + Alpha(t,m,k) * B(k,k1) * Gm(t+1,m1,k1) / NC(k1);
                end
            end
        end
    end
end

