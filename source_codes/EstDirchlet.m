function aNew = EstDirchlet(a,nx)
% Estimating single Dirichlet distribution by Newton method.
% Input:
%     a:  Initial parameter for the Dirichlet distr. (Column vector)
%     nx: Mean of log samples (Column vector).
% Output:
%     aNew: a column vector corresponding to Dirichlet parameters.

% % ---- Debug ----
% N = 5;      % Dimension. 
% T = 100;    % Number of Samples.
% Samples = rand(N,T);
% for t=1:T;Samples(:,t)=Samples(:,t)/sum(Samples(:,t));end
% nx = zeros(N,1);
% for t=1:T;nx = nx + log(Samples(:,t))/T;end
% sum(exp(nx))
% a = rand(N,1);
%------ Running's data A ----
% T = 5;
% N = 3;
% Samples =[...
%     0.0001	0.0002  0.0099  0.0096  0.0089;...
%     0.0099  0.0098  0.0001  0.0004  0.0001;...
%     0.99    0.99    0.99    0.99    0.991
%     ];
% nx = zeros(N,1);
% for t=1:T;nx = nx + log(Samples(:,t))/T;end
% a = rand(N,1);
% -------------------------
% % ---------------

N = size(nx,1);
aNew = a;
nIter = 0;
s = 10000;
obj = [];
% zrs =[norm(psi(sum(a))*repmat(1,[N,1]) - psi(a) + nx)];
while s > 1e-10
    nIter = nIter +1;
    g = psi(0,sum(aNew)) - psi(0,aNew) + nx ;
    lmd = psi(1,aNew);
    z = psi(1,sum(aNew));
    h = (g + sum(g./lmd)/(1/z - sum(1./lmd)) )./lmd;
    aNew = aNew + h;

    if find(aNew<0)
        fprintf('# Warning...Reset the parameters...\n');
        aNew = max(aNew,0.001);
        s = 10000;
    else
        s = g'*h;
    end
    %     obj = [obj log(gamma(sum(aNew)))  - sum(log(gamma(aNew))) + sum(aNew.*nx)];
    %     zrs = [zrs norm(psi(sum(aNew))*repmat(1,[N,1]) - psi(aNew) + nx)];
    
    if nIter >= 30
        fprintf('# Force to stop after 30 iterations\n');
        break;
    end
end

% obj1 = log(gamma(sum(a)))  - sum(log(gamma(a))) + sum(a.*nx);
% obj2 = log(gamma(sum(aNew)))  - sum(log(gamma(aNew))) + sum(aNew.*nx);
% obj = [obj1 obj];
% subplot(2,1,1);plot(obj)
% subplot(2,1,2);plot(zrs)
% if obj2 < obj1
%     fprintf('# Warning: Objective function decreases....\n');
% end
% obj2-obj1
% fprintf('# Number of Iters: %d\n', nIter);
