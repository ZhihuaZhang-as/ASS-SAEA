
% Create apwn sever and return handle
global aspen
aspen = actxserver('Apwn.Document.24.0'); 
% Get attributes of folder (Necessary to establish the location of the simulation)
% *.bkp file and the progam are placed in the same folder.
[stat,mess]=fileattrib; 
% Aspen Plus Simulation Name
Simulation_Name = 'Aspen_Plus_Model_for_Entrained_Flow_Coal_Gasifier';
% Open an Aspen Plus archive file and connect to the simulation engine.
aspen.invoke('InitFromArchive2',[mess.Name '\' Simulation_Name '.bkp']); %

aspen.Visible=1;           % Aspen is visible
aspen.SuppressDialogs = 1; % Suppress windows dialogs  
