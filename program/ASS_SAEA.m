function [ gbest,Record,P_Data_bestfitness] =ASS_SAEA( Data,BU,BD,problem_name,P_Data_bestfitness,fina)

% -----------------------------------------------------------------
% Important Note: This code needs intralled SURROGATE TOOLBOX(https://sites.google.com/site/srgtstoolbox/)
% -----------------------------------------------------------------
% Input:
% Data                - Data with decision variables and Objective Values
% problem_name        - Name of problem
% BU                  - Upper boundary of decision variables
% BD                  - Lower boundary of decision variables
% P_Data_bestfitness  - current best fitness
% fina                - maximum sampling multiple

% Output: 
% Record              - Logbook of evaluated solutions
% gbest               - Predicted Global best solution
% P_Data_bestfitness  - Updated current best fitness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off all;
c=size(Data,2)-1;  % NO. of decision variables
Record=[];         % Record added data 
RndData_L=[];
dis_L=[];
Class_point_L = [];
Class_p_id_L=[];
id_bound_cell_L =[];
top_bound=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StepNext=1;  %  Initializing sampling
while size(Data,1)<fina*c  % Maximum  Function Evaluations
    StepCur=StepNext;
    switch StepCur  % Current sampling
        case 1   % the most uncertainty solution
            [gbest] =PSO_U_CVV_MC(Data,BU,BD);
        case 2   % Global best solution 
            [gbest,id_M] =PSO_Gsel(Data,BU,BD);       
        case 3  % Local best solution 
            [P,ID_data]=sortrows(Data,c+1);  
            top20=floor(size(Data,1)*0.2); 
            DataT_S=P(1:top20,1:c+1);       % Top20 samples
            bu=max(DataT_S(:,1:c),[],1);    % Upper boundary
            bd=min(DataT_S(:,1:c),[],1);    % Lower boundary 
           [gbest,id_M] =PSO_Lsel(DataT_S,bu,bd);
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      num_cur=size(Data,1); % NO. of the samples
      dis=sum((abs(Data(:,1:c)-ones(num_cur,1)*gbest(:,1:c))).^2,2).^(1/2);  % Distance between samples
      [d_cur,index]=min(dis);  
      gbest(:,c+1)=compute_objectives(gbest(:,1:c),c,problem_name);   %  Function Evaluations
      best_step=[gbest,StepCur];
      Record=[Record;best_step];
       switch StepCur  
        case 1      
            StepNext=2;  % Turn to global sampling
        case 2  
            if min(Data(:,c+1))> min(gbest(:,c+1))  % If improved, keeping global sampling  
            StepNext=2;  
            else   
            StepNext=3;   % else, turning to local sampling
            end
        case 3   
            if min(Data(:,c+1))> min(gbest(:,c+1)) % If improved, keeping local sampling 
            StepNext=3;   
            else
            StepNext=1;  % else, turning to the uncertainty sampling
            end
      end
       Data=[Data;gbest(:,1:c+1)];   % Add bset solution into Data
       StepPre=StepCur;      
       P_Data_size = size(Data,1);
       P_Data_bestfitness(P_Data_size) = min(Data(:,c+1));  % Current best function values 
    end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A,I]=min(Data(:,c+1));  
gbest=Data(I,:);   

end

