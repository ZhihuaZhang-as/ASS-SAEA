function [ POP,gbest,lbest] = Update_best(POP,gbest,lbest)

% Input:
% POP           - Population of decision variables
% gbest         - Global best
% lbest         - Local best
%
% Output: 
% POP           - Updated population 
% gbest         - Updated global best
% lbest         - Updated local best
%
%------------------------------------------------------------------------
[best,Ib]=min(POP(:,end));   % min objective function
if best<=gbest(end)     
    gbest=POP(Ib,:);   
end
I=find(POP(:,end)<=lbest(:,end));  
lbest(I,:)=POP(I,:);    % optimum  in local
end