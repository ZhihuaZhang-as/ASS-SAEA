function [ gbest,id_M] =PSO_Lsel(DataL,bu,bd)

% -----------------------------------------------------------------
% Important Note: This code needs installed SURROGATE TOOLBOX(https://sites.google.com/site/srgtstoolbox/)
% -----------------------------------------------------------------
% Input:

% DataL   - Data with decision ariables and Objective Values in local region
% bu      - Upper bound of decision variables
% bd      - Lower bound of decision variables

% Output: 
%
% gbest   - Predicted Global best solution
% id_M    - ID of the selected model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=size(DataL,2)-1;      % Number of decision ariables
num_p=size(DataL,1);    % NO. of the samples in local 
dis=[]; 
Class_index=0;     % ID of sample cluster
d_p=0;  %  minimum distance
Class_point =zeros(1,c+1);  % Random samples
id=0;
id_f=0;
yhat=[];
x_best=[];
fval=[];
yhat_K=[];yhat_R=[];yhat_P=[];yhat_KRP=[];yhat_KR=[];yhat_KP=[];yhat_PR=[];
f_K=[];rank_K=[];f_R=[];rank_R=[];f_P=[];rank_P=[];f_KRP=[];rank_KRP=[];
f_KR=[];rank_KR=[];f_KP=[];rank_KP=[];f_PR=[];rank_PR=[];
rank_same=[]; yhat_mean=[];yhat_mean_s=[];
multi_x_best=[];  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Cross validation
 num_model=4;  % Total types of the models
 CVErr_model=zeros(num_p,num_model); % CV errors
 RMSE_CV=zeros(1,num_model);
for i=1:num_p
    indices(i) = i; 
end
  for i = 1:num_p 
    test = (indices == i);
    train = ~test;
    trainData =  DataL(train, :); 
    Data_x=trainData(:,1:c);   % Input of training samples
    Data_y=trainData(:,c+1);   % Output of training samples
    testData = DataL(test, :);
    x_test=testData(:,1:c);   % Input of testing samples
    y_test=testData(:,c+1);   % Output of testing samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit surrogates 
    % KRG
    srgtOPTKRG_2d  = srgtsKRGSetOptions(Data_x, Data_y);
    srgtSRGTKRG_2d = srgtsKRGFit(srgtOPTKRG_2d); 
	[yhat_K PredVar_K] = srgtsKRGPredictor(x_test, srgtSRGTKRG_2d); 
    [PRESSRMS_KRG_2d, eXV_KRG_2d] = srgtsCrossValidation(srgtOPTKRG_2d);
    % Model error
    e_y= abs((yhat_K-y_test)./y_test); 
    CVErr_model(i,1)=e_y^2;
    % Polynomial response surface
    srgtOPTPRS_2d  = srgtsPRSSetOptions(Data_x, Data_y);
    srgtSRGTPRS_2d = srgtsPRSFit(srgtOPTPRS_2d);
	[yhat_P PredVar_P] = srgtsPRSPredictor(x_test, Data_x, srgtSRGTPRS_2d); 
    [PRESSRMS_PRS_2d, eXV_PRS_2d] = srgtsCrossValidation(srgtOPTPRS_2d);
    % Model error
    e_y= abs((yhat_P-y_test)./y_test); 
    CVErr_model(i,2)=e_y^2;
    % Radial basis function
    srgtOPTRBF_2d  = srgtsRBFSetOptions(Data_x, Data_y);
    srgtSRGTRBF_2d = srgtsRBFFit(srgtOPTRBF_2d);
	yhat_R    = srgtsRBFEvaluate(x_test, srgtSRGTRBF_2d); 
    [PRESSRMS_RBF_2d, eXV_RBF_2d] = srgtsCrossValidation(srgtOPTRBF_2d);
    % Model error
    e_y= abs((yhat_R-y_test)./y_test); 
    CVErr_model(i,3)=e_y^2;
   %  KRG+RBF+PRS  computing weights
    eXVMatrix = [eXV_KRG_2d eXV_RBF_2d eXV_PRS_2d];
    CMatrix   = srgtsWASComputeCMatrix(Data_x, eXVMatrix);
    srgtsOPTs   = {srgtOPTKRG_2d  srgtOPTRBF_2d  srgtOPTPRS_2d};
    srgtsSRGTs  = {srgtSRGTKRG_2d srgtSRGTRBF_2d srgtSRGTPRS_2d};
    WAS_Model   = 'OWSdiag';
    WAS_Options = CMatrix;
    srgtOPTWAS  = srgtsWASSetOptions(srgtsOPTs, srgtsSRGTs, WAS_Model, WAS_Options);
    srgtSRGTWAS = srgtsWASFit(srgtOPTWAS); 
    w=srgtSRGTWAS.WAS_Weights;  
	[yhat_KRP PredVar_KRP] = srgtsWASPredictor(x_test, srgtSRGTWAS);   
    % Model error
    e_y= abs((yhat_KRP-y_test)./y_test); 
    CVErr_model(i,4)=e_y^2;
	end
   % RMSE_CV
   RMSE_CV=(sum(CVErr_model)*1/num_p).^(1/2);
   [RMSE_f,id_M]=min(RMSE_CV);
   switch id_M  % Selected model
   case 1 % Kriging
   srgtOPTKRG_2d  = srgtsKRGSetOptions(DataL(:,1:c), DataL(:,c+1));
   srgtSRGTKRG_2d = srgtsKRGFit(srgtOPTKRG_2d);
   % PSO 
   n=100;
   POP = initialize_sample(n,c,bu,bd);
   obj_pop = srgtsKRGPredictor(POP(:,1:c), srgtSRGTKRG_2d);
   POP=[POP,obj_pop];
   lbest=POP;
   [best,Ib]=min(POP(:,c+1));
   gbest=POP(Ib,:);
   v=(2*rand(n,c)-1).*(ones(n,1)*(bu-bd))*0.5;
   g=0;
   gmax=100;
   B=gbest;
while g<gmax  
    [ POP,v ] = fly(POP,bu,bd,gbest,lbest,v,g/gmax);
    obj_pop = srgtsKRGPredictor(POP(:,1:c), srgtSRGTKRG_2d);
    POP(:,c+1)=obj_pop;
    [ POP,gbest,lbest] = Update_best(POP,gbest,lbest);
    best=gbest(end);
    if best<B(end,end)
        B=[B;gbest];
    end
    g=g+1;
end
   % PRS
   case 2 
   srgtOPTPRS_2d  = srgtsPRSSetOptions(DataL(:,1:c), DataL(:,c+1));
   srgtSRGTPRS_2d = srgtsPRSFit(srgtOPTPRS_2d);
  % PSO 
  n=100;
  POP = initialize_sample(n,c,bu,bd);
  obj_pop = srgtsPRSPredictor(POP(:,1:c), DataL(:,1:c),srgtSRGTPRS_2d);
  POP=[POP,obj_pop];
  lbest=POP;
  [best,Ib]=min(POP(:,c+1));
  gbest=POP(Ib,:);
  v=(2*rand(n,c)-1).*(ones(n,1)*(bu-bd))*0.5;
  g=0;
  gmax=100;
  B=gbest;
while g<gmax  
    [ POP,v ] = fly(POP,bu,bd,gbest,lbest,v,g/gmax);

    obj_pop = srgtsPRSPredictor(POP(:,1:c), DataL(:,1:c),srgtSRGTPRS_2d);
    POP(:,c+1)=obj_pop;
    [ POP,gbest,lbest] = Update_best(POP,gbest,lbest);
    best=gbest(end);
    if best<B(end,end)
        B=[B;gbest];
    end
    g=g+1;
end
  % RBF
   case 3
  srgtOPTRBF_2d  = srgtsRBFSetOptions(DataL(:,1:c), DataL(:,c+1));
  srgtSRGTRBF_2d = srgtsRBFFit(srgtOPTRBF_2d);
  % PSO 
  n=100;
  POP = initialize_sample(n,c,bu,bd);
  obj_pop = srgtsRBFEvaluate(POP(:,1:c), srgtSRGTRBF_2d);
  POP=[POP,obj_pop];
  lbest=POP;
  [best,Ib]=min(POP(:,c+1));
  gbest=POP(Ib,:);
  v=(2*rand(n,c)-1).*(ones(n,1)*(bu-bd))*0.5;
  g=0;
  gmax=100;
  B=gbest;
while g<gmax  
    [ POP,v ] = fly(POP,bu,bd,gbest,lbest,v,g/gmax);

    obj_pop = srgtsRBFEvaluate(POP(:,1:c), srgtSRGTRBF_2d);
    POP(:,c+1)=obj_pop;
    [ POP,gbest,lbest] = Update_best(POP,gbest,lbest);
    best=gbest(end);
    if best<B(end,end)
        B=[B;gbest];
    end
    g=g+1;
end
   % Ensembel 
   case 4
    % KRG
    srgtOPTKRG_2d  = srgtsKRGSetOptions(DataL(:,1:c), DataL(:,c+1));
    srgtSRGTKRG_2d = srgtsKRGFit(srgtOPTKRG_2d); 
    [PRESSRMS_KRG_2d, eXV_KRG_2d] = srgtsCrossValidation(srgtOPTKRG_2d);
   
    % polynomial response surface
    srgtOPTPRS_2d  = srgtsPRSSetOptions(DataL(:,1:c), DataL(:,c+1));
    srgtSRGTPRS_2d = srgtsPRSFit(srgtOPTPRS_2d);
    [PRESSRMS_PRS_2d, eXV_PRS_2d] = srgtsCrossValidation(srgtOPTPRS_2d);
	
    % radial basis function
    srgtOPTRBF_2d  = srgtsRBFSetOptions(DataL(:,1:c), DataL(:,c+1));
    srgtSRGTRBF_2d = srgtsRBFFit(srgtOPTRBF_2d);
    [PRESSRMS_RBF_2d, eXV_RBF_2d] = srgtsCrossValidation(srgtOPTRBF_2d);
   
    eXVMatrix = [eXV_KRG_2d eXV_RBF_2d eXV_PRS_2d];
    CMatrix   = srgtsWASComputeCMatrix(DataL(:,1:c), eXVMatrix);
    srgtsOPTs   = {srgtOPTKRG_2d  srgtOPTRBF_2d  srgtOPTPRS_2d};
    srgtsSRGTs  = {srgtSRGTKRG_2d srgtSRGTRBF_2d srgtSRGTPRS_2d};
    WAS_Model   = 'OWSdiag';
    WAS_Options = CMatrix;
    srgtOPTWAS  = srgtsWASSetOptions(srgtsOPTs, srgtsSRGTs, WAS_Model, WAS_Options);
    srgtSRGTWAS = srgtsWASFit(srgtOPTWAS);   
    w=srgtSRGTWAS.WAS_Weights; 
	% PSO 
   n=100;
   POP = initialize_sample(n,c,bu,bd);
   obj_pop = srgtsWASPredictor(POP(:,1:c), srgtSRGTWAS);
   POP=[POP,obj_pop];
   lbest=POP;
   [best,Ib]=min(POP(:,c+1));
   gbest=POP(Ib,:);
   v=(2*rand(n,c)-1).*(ones(n,1)*(bu-bd))*0.5;
   g=0;
   gmax=100;
   B=gbest;
while g<gmax  
    [ POP,v ] = fly(POP,bu,bd,gbest,lbest,v,g/gmax);
    obj_pop = srgtsWASPredictor(POP(:,1:c), srgtSRGTWAS);
    POP(:,c+1)=obj_pop;
    [ POP,gbest,lbest] = Update_best(POP,gbest,lbest);
    best=gbest(end);
    if best<B(end,end)
        B=[B;gbest];
    end
    g=g+1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
