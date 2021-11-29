A dynamic model management strategy based on adaptive surrogate selection (ASS) is proposed to identify an appropriate surrogate model to assist the particle swarm optimization (PSO) algorithm faced with complex problems.      
The dynamic characteristics of the strategy are revealed in two aspects: the first is adaptive selection of the appropriate surrogate model, which can better describe the landscape of the problem in the current search space, and the second is dynamic switching between global and local modeling.                  
The most uncertain sample was also considered for searching the unexplored region and escaping from the local optimum.                 
The algorithm is verified on 10 benchmark functions and the coal gasification process.                  
The files describes the ASS-SAEA and how to communicate between MATLAB and Aspen Plus.                        

ASS_SAEA.m  
% Adopt surrogate-assisted evolutionary algorithm based on adaptive surrogate selection.         
PSO_U_CVV_MC.m  
% Find the most uncertainty solution by cross validation and Voronoi diagram.       
PSO_Lsel.m            
% Find the best solution in local region by PSO.        
PSO_Gsel.m        
% Find the best solution in global region by PSO.        
initialize_sample.m              
% Initialize samples by LHS.         
fly.m           
% Updated population in PSO algorithm.         
 Update_best.m       
% Update the best solution in PSO algorithm.       
compute_objectives.m             
% Compute the objective function of the problem.        
Create_apwn_sever.m                
% Create apwn sever and return handle.       
aspen_evaluate_objectives.m              
% Output the results of the simulation model, according to input parameter values transferred by the algorithm.     
