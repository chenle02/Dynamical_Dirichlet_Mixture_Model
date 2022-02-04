function Ind = HardClusterGamma(Gm)
% Calculate the hard clustering results from the smoothed posteriors Gamma.
% Input: 
% 	Gm: the smoothed posteriors Gamma. 1-by-D cell.
%       Gm{d} is a Td-by-M-by-K arrary: 
%         D is sequence number
%         M is the number of hidden states for Dirichlet mixture
%         K is the number of hidden states for h
%         Td is sample number of the d^th sequence 
% Output:
%   Ind: 1-by-D cell.
%       Ind{d} is a 2-by-Td matrix. The t^th sample's is generated 
%       by (h=Ind{d}(1,t),m=Ind{d}(2,t))^th Dirichlet.

Ind={};
for d = 1:length(Gm)
    Ind{d} = [];
    for t=1:size(Gm{d},1)
        k = argmax(permute(sum(Gm{d}(t,:,:),2),[3 2 1]));
        m = argmax(Gm{d}(t,:,k));
        Ind{d} = [Ind{d}, [k;m]];
    end
end

