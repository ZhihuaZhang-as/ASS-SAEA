function [ gbest] =PSO_U_CVV_MC(Data,bu,bd )

% -----------------------------------------------------------------
% Important Note: This code needs installed SURROGATE TOOLBOX(https://sites.google.com/site/srgtstoolbox/)
% -----------------------------------------------------------------
% Input:
% Data          - Data with decision ariables and objective Values 
% bu            - Upper bound of decision variables
% bd            - Lower bound of decision variables
%
% Output: 
% gbest         - Predicted global best solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=size(Data,2)-1;  %  NO. of decision variables 
x=Data(:,1:c);   % Samples
y=Data(:,c+1);   % Function Value
num_p=size(Data,1);    % NO. of samples
P_var =[];       % Variance
Class_record_p=[]; 
Class_index=0;     % ID of sample cluster
Class_point =zeros(1,c+1);  % random samples
Class_p_id=[];
bu_cell=zeros(1,c); % Upper bound of Cell
bd_cell=zeros(1,c); % Lower bound of Cell
bd_region=[];
bu_region=[];
bound_id=[];
id_bound_cell=[];  
indices=zeros(num_p,1); 
dis_region=[]; 
delt_e_sqrt=[]; 
dis=[]; % distance
test=[];
train=[];
trainData=[];
testData=[];
x_train=[];y_train=[];
x_test=[];y_test=[];
yhat=[]; PredVar=[];y_e =[];
POP=[];
PredVar_f=[];
f1=[];
MC_Class_id=[];
obj=[];
e_y= [];
dis_region_norm=[];
delt_e_sqrt_norm=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross validation
for i=1:num_p
    indices(i) = i;  
end
  for i = 1:num_p 
    test = (indices == i);
    train = ~test;
    trainData =  Data(train, :); 
    x_train=trainData(:,1:c);   % Input of training samples
    y_train=trainData(:,c+1);   % Output of training samples
    testData = Data(test, :);
    x_test=testData(:,1:c);   % Input of testing samples
    y_test=testData(:,c+1);   % Output of testing samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit surrogates
    % Kriging
    srgtOPTKRG  = srgtsKRGSetOptions(x_train, y_train);
    srgtSRGTKRG = srgtsKRGFit(srgtOPTKRG);
    [yhat PredVar] = srgtsKRGPredictor(x_test, srgtSRGTKRG);
    P_var=[P_var;x_test,PredVar,i];  
  end
   [Var_max,id_max] = max(P_var(:,c+1));   % Maximum variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   srgtOPT_all  = srgtsKRGSetOptions(x,y);
   srgtSRGT_all = srgtsKRGFit(srgtOPT_all);
% Voronoi + monte carlo sampling
w=100; 
Num_sampling = num_p*c*w;   % the number of rand samplings
RndData = zeros(Num_sampling,c); % MC samplings
for i=1:c
    RndData(:,i) =unifrnd(bd(i),bu(i),Num_sampling,1);
end
% Random points cluster
Class_record_p=zeros(Num_sampling,c+1);  
    dis=[];
for i=1:Num_sampling
    dis=sum((abs(Data(:,1:c)-ones(num_p,1)*RndData(i,1:c))).^2,2).^(1/2); % Distance function
    [d_p,Class_index]=min(dis);  
    Class_point = [RndData(i,1:c),Class_index];  
   Class_record_p(i,:) = Class_point; 
end

MC_Class_id = find(Class_record_p(:,c+1) ==id_max );  % Points in maxmum variance 
[MC_Class_num,MC_Class_num_col]= size(MC_Class_id); 
[f1 PredVar_f] = srgtsKRGPredictor(Class_record_p(MC_Class_id,1:c), srgtSRGT_all); 
dis_region=sum((abs(ones(MC_Class_num,1)*Data(id_max,1:c)-Class_record_p(MC_Class_id,1:c))).^2,2); % £¨x-xi)^2 
delt_e_sqrt=(abs(f1-ones(MC_Class_num,1)*Data(id_max,c+1))).^2;  % (f(x)-f(x*))^2
dis_min=min(dis_region);
dis_max=max(dis_region);
delt_e_sqrt_min=min(delt_e_sqrt);
delt_e_sqrt_max=max(delt_e_sqrt);
dis_region_norm=(dis_region-ones(MC_Class_num,1)*dis_min)/(dis_max-dis_min); % Normalization
delt_e_sqrt_norm=(delt_e_sqrt-ones(MC_Class_num,1)*delt_e_sqrt_min)/(delt_e_sqrt_max-delt_e_sqrt_min);
obj=dis_region_norm.*delt_e_sqrt_norm;  % the uncertainty of objective function

[best,Ib]=max(obj);
gbest=Class_record_p(MC_Class_id(Ib),:);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end