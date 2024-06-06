%%%%%%%%%%%%%%%%%%%%%%%
%
% Code to calculate the sigantures of metastability described in 
% the NRN paper 'MAetastability demystified
% Doi:
%
% Fran Hancock, June 2024
%
% This code is a collection of snipits taken from publically available
% code. References are provided.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%

clc

p = path;
restoredefaultpath % need to do this to clear my extensive path

wmi = pwd;
addpath(genpath(wmi));

%%%%%% read in the defaults for this project

param = loadParameters; 

DRV_DATA = param.DRV_DATA;
MAT_DATA = param.MAT_DATA;
NSUBS = param.NSUBS;


N_areas = param.N;
Rmax = param.Rmax;

TR = param.TR;

Tmax = param.Tend; % V1_all and Time_sessions are reduced to actual TRs


%%   
% LEADING EIGENVECTOR DYNAMICS ANALYSIS
%
% Written by Joana Cabral, May 2016 (joanacabral@med.uminho.pt)
% Modified by Joana Cabral, November 2020
% Modified by Fran Hancock, November 2021
%
% Step 1 -  read in the parcellated files
%           Bandpass filter the fMRI images
%           Apply the Hilbert transform to obtain the analytical signal
%           Calculate the Kuramoto order parameter and the standard deviation
%           Calculate the pairwise phase difference
%           Obtain the eigenvector and eigenvalue for each TR
%           Save these in LEiDA_EigenVectors in V1_all, val1_all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eigs_file = 'LEiDA_EigenVectors';
Sigs_file = 'ALL_Sigs';

for run = 1:Rmax 
	
    status = mkdir(['RUN' num2str(run)]);
    
    Data_info = dir([MAT_DATA num2str(run) '/AAL*Sub-*']);
    % Define total number of scans that will be read
    total_scans=size(Data_info,1);

    if total_scans ~= 0
    
        % Create empty variables to save patterns 
        
        V1_all = zeros(total_scans*(Tmax-2),N_areas); % Alerady in shape for the kmeans
        val1_all = zeros(total_scans*(Tmax-2),1);

        %val1_all=zeros(total_scans*(Tmax-2)); % Alerady in shape for the kmeans
        Time_sessions=zeros(1,total_scans*(Tmax-2)); % to know the scan number 
        t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*(Tmax-2)) 
        
        for sub=1:total_scans
            
             disp(Data_info(sub).name)
            
            % load BOLD signals in each scan
            signal = struct2array(load([MAT_DATA num2str(run) '/' Data_info(sub).name]));
            
            %%%%%%%%%%%%%%% Band pass filter before getting the analtical signal %%%%%%%%%%%%%%%%%
            
            high_pass=0.01;     % Bands from Glerean but with fast TR
            low_pass=0.08;
        
            for n=1:size(signal,1)
                signal(n,:)=bandpass(signal(n,:),[high_pass low_pass],1/TR);
            end
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Get the fMRI phase using the Hilbert transform

            for seed=1:N_areas
                ts=signal(seed,:);
                signal(seed,:) = angle(hilbert(ts));  % angle returns the instantaneous phase from the analytical signal 
            end
            
            % now calculate the KOP
            GLOBAL_KOP(sub,:) = abs(sum(exp(1i*signal)))/N_areas;% abs(mean(sum of e itheta(t)))
            GLOBAL_META(sub,:) = std(GLOBAL_KOP(sub,:),1);

            GLOBAL_IGNITE(sub,:) = get_global_ignite(signal, N_areas);

            for t=2:size(signal,2)-1 % skip the first and last values
                
                [v1,val1]=eigs(cos(signal(:,t)-signal(:,t)'),1);  % get the 1st eigenvector and eigenvalue
                
                if sum(v1)>0  % ensure the 1st is negative - this is just a convention
                    v1=-v1;
                end
                
                t_all=t_all+1;
                V1_all(t_all,:)=v1;
                val1_all(t_all)=val1;
                Time_sessions(t_all)=sub;
            end
        end
        
        
        % Reduce size in case some scans have less TRs than Tmax
        V1_all(t_all+1:end,:)=[];
        val1_all(t_all+1:end,:)=[];
        Time_sessions(:,t_all+1:end)=[];

        val1_all = val1_all';

  
    else
         disp([' Could not find the parcellated files for ' MAT_DATA num2str(run) '/AAL*sub-*'])
    end
    
    save([DRV_DATA num2str(run) '/' eigs_file],'V1_all','val1_all','Time_sessions')
    save([DRV_DATA num2str(run) '/' Sigs_file],'GLOBAL_KOP','GLOBAL_META','GLOBAL_IGNITE')  

    disp([' Eigenvectors and eigenvalues saved successfully as ' DRV_DATA num2str(run) '/' eigs_file])

    clear val1_all

end


%%
% now calculate 1) the mean VAR signature, and the spectral radius
% signature
% 

for run = 1:Rmax

    load([DRV_DATA num2str(run) '/LEiDA_EigenVectors.mat']); % V1_all, val1_all; Time_sessions
    
    

    for sub = 1:NSUBS(run)

        sub_idx = Time_sessions == sub; % pick out the entries for this subject
        GLOBAL_VAR_MEAN(sub,:) = mean(var(V1_all(sub_idx,:),1));
        GLOBAL_SR(sub,:) = std(val1_all(:,sub_idx));

    end

    save([DRV_DATA num2str(run) '/' Sigs_file],'GLOBAL_VAR_MEAN','GLOBAL_SR', '-append')

end


