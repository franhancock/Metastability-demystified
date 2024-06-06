classdef loadParameters < matlab.mixin.Copyable
%% loadParameters.m     
%
% Contains all the parameters of the imaging data for metric calculations. It is
% necessary to make an instance of the class before you can use the 
% parameters.
%
% Example:
% >> param = utils.loadParameters;
% 
%
% Some important notes:
% 1. All the parameters can be changed by overwriting the existing instance 
%
% 2. The methods section dynamically calculates/refreshes the dependent
%    parameters when other independent parameters are changed.
%
% Original: James Pang, QIMR Berghofer, 2020
%
% Adapted for COBRE ANALYTICAL COMPLEXITY
% Fran Hancock 2023
%
%%
    properties 
    % =====================================================================
    %                IMAGING AND FOLDER PARAMETERS
    % ===================================================================== 
        
        TR              = 2;                        % connectivity matrix
        Tstart      	= 1;
        Tend       		= 150;
        N               = 116;                      % number of nodes 
                                                    
        NSUBS           =[20 20];  

        STUB_FILE       = 'AAL_Sub-00'              % parcellated file name stub
       
        MAT_DATA        = 'COBRE_MAT_R';            % parcellated data Set to 
        DRV_DATA        = 'RUN';                    % derivatives    

        Rmax            = 2;                        % number of runs/conditions
        parcellation    = 'AAL116';                 % Parcellation

        cond            = ['CONT';'SCHZ'];
        proj            = "COBRE";

         alpha           = 0.85; 

    end
end
       