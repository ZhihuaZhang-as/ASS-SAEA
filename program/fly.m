function [ POP,v ] = fly(POP,bu,bd,gbest,lbest,v,theta)

% Input:
% POP           - Population of decision variables
% gbest         - Global best
% lbest         - Local best
% bu            - Upper boundary of decision variables
% bd            - Lower boundary of decision variables
% v             - Velocity
% theta         - Generation parameter
%               
% Output: 
% POP           - Updated population 
% v             - Updated velocity
%
%------------------------------------------------------------------------
c1=1.49445;
c2=1.49445;
N=size(POP,1);  % NO. of population
n=size(bu,2);   % NO. of decision variables
w=0.9-theta*0.5;  % Updated inertia weight
v=w*v+c1*rand(N,n).*(lbest(:,1:n)-POP(:,1:n))+c2*rand(N,n).*(ones(N,1)*gbest(1:n)-POP(:,1:n));  %  Updated velocity
tvmax=0.5*ones(N,1)*(bu-bd);   % Upper boundary of velocity
I=find(v>tvmax);   
v(I)=tvmax(I);    
I=find(v<(0-tvmax)); 
v(I)=0-tvmax(I);
tPOP=POP(:,1:n)+v;  % Updated position
Tbd=ones(N,1)*bd;   
Tbu=ones(N,1)*bu;   
I=find(tPOP>Tbu);  
 tPOP(I)=Tbu(I)-(tPOP(I)-Tbu(I));   %  

I=find(tPOP<Tbd);
 tPOP(I)=Tbd(I)+(Tbd(I)-tPOP(I));  %    

POP(:,1:n)=tPOP;
end