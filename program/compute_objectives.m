function [ obj] =compute_objectives(POP,c,problem_name)

% Input:
% problem_name  - Benchmark problem
% c             -No. of Decision Variables
% POP           -Population of Decision Variables
%
% Output: 
% obj           - Calculated objective value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj=[];
n=size(POP,1);   % NO. of Population
obj_m1=0;
obj_m2=0;
x=[];
y=[];
switch problem_name   
    case 'Rastrigin'   % minf=0 x=[0 0...0];  Bounds[-5,5]   
        obj=10*c+sum(POP(:,1:c).^2-10*cos(2*pi*POP(:,1:c)),2);
    case 'Rosenbrock'  % min=0,x=[1 1...1];  Bounds[-2.048£¬2.048] 
        P=100*(POP(:,2:c)-POP(:,1:c-1).^2).^2;
        obj=sum(P,2)+sum((POP(:,1:c-1)-1).^2,2); 
    case 'Griewank'  % min=0; x=[0 0...0];Bounds ¡Ê [-600£¬600] 
        P=[1:c].^0.5;
        P=ones(size(POP,1),1)*P;
        obj=sum(POP(:,1:c).^2,2)/4000+1-prod(cos(POP(:,1:c)./P),2);
    case 'Ellipsoid'  % min=0; x=[0 0...0];Bounds[-5.12,5.12]  
        P=[1:c];
        P=ones(size(POP,1),1)*P;
        obj=sum(P.*(POP(:,1:c).^2),2); 	
    case 'Schwefel'      %Bounds[-500,500]  
          obj=418.9828872724338*c-sum(POP(:,1:c).*sin((abs(POP(:,1:c))).^0.5),2);   
end

end

