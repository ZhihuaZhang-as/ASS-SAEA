function [ init_sample ] = initialize_sample(n,c,bu,bd)

% Input:
% bu       - Upper bound
% bd       - Lower bound
% c        - No. of decision variables
% n        - Samples scale
%
% Output: 
% init_sample    -Initial samples

init_sample=lhsdesign(n,c).*(ones(n,1)*(bu-bd))+ones(n,1)*bd;    

end