%%
% Now calculate the standard deviation of intrinsic iggnition events
%
% Code from Gustavo Deco
% Eddited by Sonsoles Alonso Martinez
% alonsomartinez.sonsoles@gmail.com
% Last edited Feb 2019
%
% Adapted for the Metastability demystified paper
% Fran Hancock
% fran.hancock@kcl.ac.uk
%
% June 2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function global_ignite = get_global_ignite(sub_signal, N_areas)

Isubdiag = find(tril(ones(N_areas),-1));

nTRs = 3; % nTRs size window

        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Account for less than 152 timepoints
            Tact=size(sub_signal,2);
            T = 1:Tact; 
            Tmax_sub = Tact;

            events=zeros(N_areas,Tmax_sub); 

            for seed=1:N_areas
                ts = sub_signal(seed,:); 
                tise = detrend(demean(ts(T)));

                ev1 = tise>std(tise)+mean(tise);
                ev2 = [0 ev1(1:end-1)];
                events(seed,:) = (ev1-ev2)>0;     
            end
 
    
            % Get meassures
               
            % ----------------------- integration --------------------------------
            integ=zeros(1,T(end));
            for t = T
                % For each time point a matrix of co-occurring events is calculated
                events_matrix=zeros(N_areas,N_areas);
                for i = 1:N_areas 
                    for j = 1:N_areas 
                      events_matrix(i,j) = events(i,t)*events(j,t); 
                    end
                end 
                A = events_matrix; 
                A = A-eye(N_areas);     
                [~, csize] = get_components(abs(A));% number of regions within comp.
                cs = max(csize); % size of the largest component
                integ(t)=cs/N_areas;  
            end 
        
            % --------------------- ignition-driven integration -----------------
            nevents2 = zeros(N_areas,1);
            IntegStim2 = zeros(N_areas,N_areas,1);
            % save events and integration values for nTRs after the event
            for seed = 1:N_areas
                  flag = 0; % flag 0 to start counting the events
                  for t = T
                        %detect first event (nevents = matrix with 1xnode and
                        %number of events in each cell)
                        if events(seed,t) == 1 && flag == 0 
                          flag = 1;
                          nevents2(seed) = nevents2(seed)+1;  
                        end
                        %save integration value for nTRs after the first event
                        %(nodesx(nTR-1)xevents)
                        if flag > 0
                          IntegStim2(seed,flag,nevents2(seed)) = integ(t);
                          flag = flag+1;
                        end
                        %after nTRs, set flag to 0 and wait for the next event
                        %(then, integ saved for (nTRs-1) events)
                        if flag == nTRs
                          flag = 0;
                        end
                  end
            end
        
        %---------------------------------------------------------------
        %                       LOCAL MEASSURES
        %---------------------------------------------------------------
        stdevokedinteg2 = zeros(N_areas,1);
        
        for seed = 1:N_areas
          stdevokedinteg2(seed) = std(max(squeeze(IntegStim2(seed,:,1:nevents2(seed))))); 
        end
    
        %---------------------------------------------------------------
        %                       GLOBAL MEASSURES
        %---------------------------------------------------------------    
        global_ignite = mean(stdevokedinteg2);
end
