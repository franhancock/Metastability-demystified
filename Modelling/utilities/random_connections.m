function Coupling = random_connections
% 
% function to generate 32 random connections from each community of 32
% oscillators to any to any other community
%
% Replicate Shanahan (2010) model
%
% Fran Hancock, 2023
% fran.hancock@kcl.ac.uk
%

N_oscillators = 32;
N_communities = 8;
Total_oscillators = N_communities*N_oscillators;

u = 0.6; % intracommunity coupling strength
v = 0.4; % intercommunitiy coupling strength

kij = zeros(Total_oscillators, Total_oscillators);
% Create the intracommunity couplings
for communities=1:N_communities
    if communities == 1
        start_community = 1;
        kij(start_community:N_oscillators,start_community:N_oscillators)=ones(N_oscillators);     
    else
        start_community = (communities-1)*N_oscillators;
        kij(start_community+1:start_community+N_oscillators,start_community+1:start_community+N_oscillators)=ones(N_oscillators);     
    end

    Kij = kij*u;
end
Kij=Kij-diag(diag(Kij));

% now make 32 random connections from each oscilator in a community to any
% other oscillators

%random_links(:,:) = randi([0,Total_oscillators],N_oscillators,N_oscillators);

for communities=1:N_communities
    if communities == 1

        nodes=[33:256];
        
        % now assign coupling to the nodes
        for node=1:32
            % get 32 random nodes to connect to
            random_nodes=randsample(nodes,N_oscillators);
            
            for coupling=1:32
                % disp([num2str(random_nodes(coupling)) ' ' num2str(node)])
                Kij(random_nodes(coupling), node) = v;
                Kij(node,random_nodes(coupling)) = v;
            end
        end   
    else
        nodes=[1:N_oscillators*(communities-1) communities*N_oscillators+1:Total_oscillators];
     
        % now assign coupling to the nodes
        for node=N_oscillators*(communities-1)+1:communities*N_oscillators
            
            random_nodes=randsample(nodes,N_oscillators);

            for coupling=1:32
%               disp([num2str(random_nodes(coupling)) ' ' num2str(node)])
                Kij(random_nodes(coupling), node) = v;
                Kij(node,random_nodes(coupling)) = v;

%                Kij(random_nodes(node - (communities-1)*N_oscillators), node) = v;
%                Kij(node,random_nodes(node - (communities-1)*N_oscillators)) = v;
           end
        end   
    end

end
Coupling = Kij;

