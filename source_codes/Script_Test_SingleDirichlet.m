% Demo and test script for EstDirchlet.m

% --- Parameters -----
M = 1;
T =2000;
a0 = [3 2 4 5 8 10 20];
% --------------------

N = size(a0,2);
DataSet = GenDir(a0', T);
A = MomentMatchingInitDM(M,DataSet);
gm = T;
nx = zeros(N,1);
for t=1:T
    nx = nx + log(DataSet(:,t));
end
nx = nx/gm;
aNew = EstDirchlet(A',nx);

dataLikelihood_Ori = sum(dataLikelihood_DM(a0,DataSet,1))/T;
dataLikelihood_MM = sum(dataLikelihood_DM(A,DataSet,1))/T;
dataLikelihood_EST =sum(dataLikelihood_DM(aNew',DataSet,1))/T;

[a0' A' aNew;...
    0 norm(a0-A) norm(a0-aNew');...
    dataLikelihood_Ori dataLikelihood_MM dataLikelihood_EST]

% columnLabels = {'Real Parameter', 'Para. by Moment Matching','Estimated Para.'};
% rowLabels ={};
% for n=1:N
%     rowLabels = [rowLabels strcat('Dim. ',int2str(n))];
% end
% rowLabels = [rowLabels,'Euclidian dist.','Average data log likelihood'];
% matrix2latex(...
%     [a0' A' aNew;...
%     0 norm(a0-A) norm(a0-aNew');...
%     dataLikelihood_Ori dataLikelihood_MM dataLikelihood_EST],...
%     'lele.tex',...
%     'columnLabels', columnLabels,...
%     'rowLabels', rowLabels,...
%     'alignment', 'c',...
%     'format', '%-6.4f')

