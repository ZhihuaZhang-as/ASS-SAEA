 function [EFFECT_GAS_PRO,Temp, C_conv, EFFECT_GAS]  = aspen_evaluate_objectives(X, V,aspen)

% Input:
% X                 - Decision variables
% V                 - NO.of decision variables

% Output:
% EFFECT_GAS_PRO    - the objective functions of the effective gas productivity
% Temp              - Temparature
% C_conv            - Carbon conversion
% EFFECT_GAS        - the mole percent of the effective gas

% Here, the algorithm maximizes the objective function hence 
% if you would like to minimizes the function then multiply the function by negative one¡£
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_p=size(X,1);    % the number of samples
EFFECT_GAS_PRO=zeros(num_p,1);
EFFECT_GAS=zeros(num_p,1);
Temp=zeros(num_p,1);
C_conv=zeros(num_p,1);
for i= 1:num_p
% Transfer data in MATLAB into Aspen Plus
aspen.Tree.FindNode('\Data\Streams\COAL\Input\FLOW\NC\COAL').Value=X(i,1);  % Coal mass-flow 
aspen.Tree.FindNode('\Data\Streams\COAL\Input\TEMP\NC').Value=X(i,2);  % Coal temparature 
aspen.Tree.FindNode('\Data\Streams\COAL\Input\PRES\NC').Value=24;  %  Coal  pressure
aspen.Tree.FindNode('\Data\Streams\OXYGEN\Input\TOTFLOW\MIXED').Value=X(i,1)*X(i,3);  % O2 total-flow 
aspen.Tree.FindNode('\Data\Streams\OXYGEN\Input\TEMP\MIXED').Value=X(i,6);    % O2 temparature
aspen.Tree.FindNode('\Data\Streams\OXYGEN\Input\PRES\MIXED').Value=24;    %  O2 pressure  
aspen.Tree.FindNode('\Data\Streams\STEAM\Input\TOTFLOW\MIXED').Value=X(i,1)*X(i,4);   % Steam total-flow
aspen.Tree.FindNode('\Data\Streams\STEAM\Input\TEMP\MIXED').Value=X(i,5);     % Steam temparature
aspen.Tree.FindNode('\Data\Streams\STEAM\Input\PRES\MIXED').Value=24;        % Steam pressure
aspen.Tree.FindNode('\Data\Blocks\GASIFIER\Input\PRES').Value=24;          % Gasifier pressure
aspen.Tree.FindNode('\Data\Blocks\COMBUST\Input\PRES').Value=24;          %  COMBUST pressure
aspen.Tree.FindNode('\Data\Blocks\MIXER\Input\PRES').Value=24;         % MIXER pressure
aspen.Tree.FindNode('\Data\Blocks\PRESCORR\Input\PRES').Value=24;     %  PRESCORR pressure
aspen.Tree.FindNode('\Data\Blocks\SEPELEM\Input\PRES').Value=24;     %  SEPELEM pressure
aspen.Tree.FindNode('\Data\Blocks\SEPSG\Input\PRES1').Value=24;     %   SEPSGpressure


Reinit(aspen) 
aspen.Run2
% Receive the output of Aspen Plus
% Gas flow kmol/sec  
 try
 CO = str2double(aspen.Tree.FindNode('\Data\Streams\PRODUCT\Stream Results\Table\(  CO,PRODUCT)').Value);                          
 H2 = str2double(aspen.Tree.FindNode('\Data\Streams\PRODUCT\Stream Results\Table\(  H2,PRODUCT)').Value);  
 CH4= str2double(aspen.Tree.FindNode('\Data\Streams\PRODUCT\Stream Results\Table\(  CH4,PRODUCT)').Value);  
 CO2= str2double(aspen.Tree.FindNode('\Data\Streams\PRODUCT\Stream Results\Table\(  CO2,PRODUCT)').Value);  
 H2S= str2double(aspen.Tree.FindNode('\Data\Streams\PRODUCT\Stream Results\Table\(  H2S,PRODUCT)').Value);  
 N2=  str2double(aspen.Tree.FindNode('\Data\Streams\PRODUCT\Stream Results\Table\(  N2,PRODUCT)').Value);  
 C6H6=str2double(aspen.Tree.FindNode('\Data\Streams\PRODUCT\Stream Results\Table\(  C6H6,PRODUCT)').Value);  
 % Gasifier temparature
 Temp(i)=str2double(aspen.Tree.FindNode('\Data\Streams\PRODUCT\Stream Results\Table\(Temperature K,PRODUCT)').Value);      
 % Calculate carbon conversion
 COAL_H2O = aspen.Tree.FindNode('\Data\Streams\COAL\Input\ELEM\NC\COAL\PROXANAL\#0').Value; % coal water percent
 COAL_C = aspen.Tree.FindNode('\Data\Streams\COAL\Input\ELEM\NC\COAL\ULTANAL\#1').Value; % ULTANAL  carbon percent
 C_conv(i)=(CO+CH4+CO2+C6H6*6)*12*1000/(X(i,1)*(1-COAL_H2O*0.01)*COAL_C*0.01);   
 % Effective gas productivity
 EFFECT_GAS_PRO(i)= (CO+H2)*1000*22.4*0.001./(X(i,1)*0.001);  % unit m3/kg
 % Effective gas 
 EFFECT_GAS(i)= (CO+H2)/(CO+H2+CH4+CO2+H2S+N2);  
 catch  % if the calcualation is abnormal, it ignores the sample.
     EFFECT_GAS_PRO(i)=0.1;
      Temp(i)=0.1;
      C_conv(i)=0.1;
 end
 
end

% aspen.Close;
% aspen.Quit;
 end 

 