function [ gbest,id_M] =PSO_Gsel(Data,bu,bd )

% -----------------------------------------------------------------
% Important Note: This code needs installed SURROGATE TOOLBOX(https://sites.google.com/site/srgtstoolbox/)
% -----------------------------------------------------------------
% Input:
% Data          - Data with decision variables and Objective Values
% bu            - Upper boundary of decision variables
% bd            - Lower boundary of decision variables
%
% Output: 
 
% gbest         - Predicted Global best solution
% id_M          - ID of the selected model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=size(Data,2)-1; %  NO. of decision ariables
x=Data(:,1:c);   % Samples
y=Data(:,c+1);   % Function Value
num_p=size(Data,1);  % NO. of samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross validation
num_model=4; % Total types of the models
kfold=10;
indices=[];
%  K-fold cross validation
indices=crossvalind('Kfold',num_p,kfold);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVErr_model=zeros(kfold,num_model);  % Cross validation
RMSE_CV=zeros(1,num_model);

for i=1:kfold
    test =[];
    train = [];
    trainData = [];
    testData = [];
    Data_x=[];
     Data_y=[];
     x_test=[];
     y_test=[];
     test = (indices == i);
    train = ~test;
    trainData =  Data(train, :); 
    Data_x=trainData(:,1:c);  % Input of training samples
    Data_y=trainData(:,c+1);  % Output of training samples
    testData = Data(test, :);
    x_test=testData(:,1:c);  % Input of testing samples
    y_test=testData(:,c+1);   % Output of testing samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 polypart=mean(Data_y);
    % Fit surrogates 
    % KRG
    srgtOPTKRG_2d  = srgtsKRGSetOptions(Data_x, Data_y);
    srgtSRGTKRG_2d = srgtsKRGFit(srgtOPTKRG_2d); 
	[yhat_K PredVar_K] = srgtsKRGPredictor(x_test, srgtSRGTKRG_2d); 
    [PRESSRMS_KRG_2d, eXV_KRG_2d] = srgtsCrossValidation(srgtOPTKRG_2d);
    e_y= abs((yhat_K-y_test)./y_test); 
    CVErr_model(i,1)=sum(e_y.^2);

    % polynomial response surface
    srgtOPTPRS_2d  = srgtsPRSSetOptions(Data_x, Data_y);
    srgtSRGTPRS_2d = srgtsPRSFit(srgtOPTPRS_2d);
	[yhat_P PredVar_P] = srgtsPRSPredictor(x_test, Data_x, srgtSRGTPRS_2d); 
    [PRESSRMS_PRS_2d, eXV_PRS_2d] = srgtsCrossValidation(srgtOPTPRS_2d);

    e_y= abs((yhat_P-y_test)./y_test); 
    CVErr_model(i,2)=sum(e_y.^2);

    srgtOPTRBF_2d  = srgtsRBFSetOptions(Data_x, Data_y);
    srgtSRGTRBF_2d = srgtsRBFFit(srgtOPTRBF_2d);
	yhat_R    = srgtsRBFEvaluate(x_test, srgtSRGTRBF_2d); 
    [PRESSRMS_RBF_2d, eXV_RBF_2d] = srgtsCrossValidation(srgtOPTRBF_2d);

    e_y= abs((yhat_R-y_test)./y_test); 
    CVErr_model(i,3)=sum(e_y.^2);

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

    e_y= abs((yhat_KRP-y_test)./y_test); 
    CVErr_model(i,4)=sum(e_y.^2);
end
   % RMSE_CV
   RMSE_CV=(sum(CVErr_model)/num_p).^(1/2); 
   [RMSE_f,id_M]=min(RMSE_CV);

switch id_M   % Selected model
   case 1
  % Kriging
    srgtOPTKRG_2d  = srgtsKRGSetOptions(Data(:,1:c), Data(:,c+1));
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
  srgtOPTPRS_2d  = srgtsPRSSetOptions(Data(:,1:c), Data(:,c+1));
  srgtSRGTPRS_2d = srgtsPRSFit(srgtOPTPRS_2d);
  % PSO 
  n=100;
  POP = initialize_sample(n,c,bu,bd);
  obj_pop = srgtsPRSPredictor(POP(:,1:c), Data(:,1:c),srgtSRGTPRS_2d);
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
    obj_pop = srgtsPRSPredictor(POP(:,1:c), Data(:,1:c),srgtSRGTPRS_2d);
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
   srgtOPTRBF_2d  = srgtsRBFSetOptions(Data(:,1:c), Data(:,c+1));
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

    srgtOPTKRG_2d  = srgtsKRGSetOptions(Data(:,1:c), Data(:,c+1));
    srgtSRGTKRG_2d = srgtsKRGFit(srgtOPTKRG_2d); 
    [PRESSRMS_KRG_2d, eXV_KRG_2d] = srgtsCrossValidation(srgtOPTKRG_2d);
   
    % polynomial response surface
    srgtOPTPRS_2d  = srgtsPRSSetOptions(Data(:,1:c), Data(:,c+1));
    srgtSRGTPRS_2d = srgtsPRSFit(srgtOPTPRS_2d);
    [PRESSRMS_PRS_2d, eXV_PRS_2d] = srgtsCrossValidation(srgtOPTPRS_2d);
	
    srgtOPTRBF_2d  = srgtsRBFSetOptions(Data(:,1:c), Data(:,c+1));
    srgtSRGTRBF_2d = srgtsRBFFit(srgtOPTRBF_2d);
    [PRESSRMS_RBF_2d, eXV_RBF_2d] = srgtsCrossValidation(srgtOPTRBF_2d);
   
    eXVMatrix = [eXV_KRG_2d eXV_RBF_2d eXV_PRS_2d];
    CMatrix   = srgtsWASComputeCMatrix(Data(:,1:c), eXVMatrix);
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

end








