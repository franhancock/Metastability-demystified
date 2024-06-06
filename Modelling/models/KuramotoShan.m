% KuramotoShan - Network of Kuramoto Phase Oscillators 
%   Constructs a Kuramoto network with n nodes.
%       theta_i' = omega_i + SUM_j Kij*sin(theta_i-theta_j) - alpha
%   where 
%       theta is an (nx1) vector of oscillator phases (in radians),
%       omega is an (nx1) vector of natural frequencies (cycles/sec)
%       Kij is an (nxn) matrix of connection weights,
%
%
% Adapted by Fran Hancock 2023 from:
%
% Authors
%   Stewart Heitmann (2016a,2017a,2018a,2018b,2020a)
% --------------------------------------------------
%
% For the model in Shanahan 2010
%
% N = 63 - total of intra and inter couplings
% Kij - A coupling matrix with full intracoupling and random intercoupling
% u = 0.6 intracoupling strength
% v = 0.4 intercoupling strength
% Alpha = 1.47 - phase lag 
% Omega = 1 Oscillator frequencies
% Theta = random initial phases
% N_communities = 8
% N_oscillators = 32;
%
% Call:
% 
% Kij = random_connections; sys = KuramotoShan(Kij); gui = bdGUI(sys);  
%
% Copyright (C) 2016-2022 QIMR Berghofer Medical Research Institute
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in
%    the documentation and/or other materials provided with the
%    distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
% COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% generate 32 random connections from each community of 32
% oscillators to any to any other community


function sys = KuramotoShan(Kij)
    addpath(genpath('/Users/HDF/Hancock Dean Dropbox/Doug Dean/Fran/Academics/PhD/Matlab_World/bdtoolbox-2022b'))
    % Determine the number of nodes from the size of the coupling matrix
    n=size(Kij,1);
%     N = n; %63;
%     N_communities = 8;
%     N_oscillators=32;

    % Handle to our ODE function
    sys.odefun = @odefun;
    
    % ODE parameters
    sys.pardef = [ struct('name','Kij',   'value',Kij);
                   struct('name','omega', 'value',ones(n,1));
                   struct('name','Alpha', 'value',1.47);                                                                                                           
                   ];
      
    % ODE state variables
   sys.vardef = struct('name','theta', 'value',2*pi*rand(n,1), 'lim',[0 2*pi]);
    
    % Time span
    sys.tspan = [0 1000];
    sys.odeoption.MaxStep = 0.05; 
    sys.tstep = 0.05;
    

    % Relevant ODE solvers
    sys.odesolver = {@ode45,@ode23,@ode113,@odeEul};
    
    % ODE solver options
    sys.odeoption.RelTol = 1e-6;
                   
                    
    % Latex (Equations) panel
    sys.panels.bdLatexPanel.title = 'Equations'; 
    sys.panels.bdLatexPanel.latex = {
        '$\textbf{KuramotoShan}$'
        ''
        'A community coupled network of Kuramoto Oscillators'
        '{ }{ }{ } $\dot \theta_i = \omega_i - \frac{1}{N} \sum_j K_{ij} \sin(\theta_i - \theta_j - \alpha)$'
        'where'
        '{ }{ }{ } $\theta_i$ is the phase of the $i^{th}$ oscillator (radians),'
        '{ }{ }{ } $\omega_i$ is its natural oscillation frequency (cycles/sec),'
        '{ }{ }{ } $K$ is the network connectivity matrix ($n$ x $n$),'
        '{ }{ }{ } $alpha$ is is a fixed phase lag,'
        '{ }{ }{ } $i,j=1 \dots n,$'
        ['{ }{ }{ } $n{=}',num2str(n),'.$']
        ''
        'The Kuramoto order parameter ($R$) is a metric of phase synchronisation.'
        '{ }{ }{ } $R = \frac{1}{N_c} \| \sum_i \exp(\mathbf{i} \theta_i) \|$'
        'It corresponds to the radius of the centroid of the phases, as shown in'
        'the Auxiliary panel.'
        ''
        'References'
        'Kuramoto (1984) Chemical oscillations, waves and turbulence.'
        'Strogatz (2000) From Kuramoto to Crawford.'
        'Shanahan (2010) Metastable chimera states in community-structured oscillator networks.'
        'Breakspear et al (2010) Generative models of cortical oscillations.'
        'Chapter 6.2 of the Handbook for the Brain Dynamics Toolbox (Version 2022b).'
        };
    
    % Time Portrait panel
    sys.panels.bdTimePortrait.title = 'Time Portrait';
    sys.panels.bdTimePortrait.modulo = 'on';
 
    % Phase Portrait panel
    sys.panels.bdPhasePortrait.title = 'Phase Portrait';
    sys.panels.bdPhasePortrait.modulo = 'on';

    % Auxiliary panel
    %sys.panels.bdAuxiliary.title = 'Auxiliary';
    sys.panels.bdAuxiliary.auxfun = {@centroid1,@KuramotoR};


    % Solver panel
    sys.panels.bdSolverPanel.title = 'Solver';                
end

% Kuramoto ODE function where
% theta is a (1xn) vector of oscillator phases,
% Kij is either a scalar or an (nxn) matrix of connection weights,
% omega is either a scalar or (1xn) vector of oscillator frequencies.
function dtheta = odefun(t,theta,Kij,omega,Alpha)
    n = numel(theta);
    theta_i = theta * ones(1,n);                        % (nxn) matrix with same theta values in each row
    theta_j = ones(n,1) * theta';                       % (nxn) matrix with same theta values in each col
    theta_ij = theta_i - theta_j;                       % (nxn) matrix of all possible (theta_i - theta_j) combinations
    dtheta = omega + 1/n.*sum(Kij.*sin(theta_ij-Alpha),1)';   % Kuramoto Equation in vector form.

end

% Auxiliary function that plots the centroid of the oscillators
function centroid1(ax,t,sol,Kij,omega,Alpha)

    N_communities = 8;
    N_oscillators = 32;

    % Get the phases of the oscillators at time t
    theta = bdEval(sol,t);
    
    ztheta = exp(1i.*theta);
    
    % Get the centroid for each community
    for c=1:N_communities
        if c==1
            ctheta(c,:)=sum(ztheta(1:32,:));
            mu(c,:) = angle(ctheta(c,:));
        else
            ctheta(c,:)=sum(ztheta((c-1)*N_oscillators:c*N_oscillators,:));
            mu(c,:) = angle(ctheta(c,:));
        end
    end

    ctheta = exp(1i.*mu);
    % Project the phases into the complex plane.


    % Plot the centroid.
    centroidplot(ax,ctheta,ztheta);
    %text(ax,-1,-1,num2str(t,'time = %g'));
    title(ax,'centroid of oscillators'); 
end


% Utility function for plotting the centroid
function centroidplot(ax,ztheta,atheta)
    % compute the phase centroid
    centroid = mean(ztheta);
    
    % plot the unit circle
    plot(ax,exp(1i.*linspace(-pi,pi,100)), 'color',[0.75 0.75 0.75]);
    
    % plot the oscillator phases on the unit circle
    if size(ztheta,1)==8
        plot(ax,[0 ztheta(1)], 'color', 'b');
        plot(ax,ztheta(1),'o','color','k','MarkerFaceColor','b', 'MarkerSize',14);
         plot(ax,[0 ztheta(2)], 'color', 'g');
        plot(ax,ztheta(2),'o','color','k','MarkerFaceColor','g', 'MarkerSize',14);
         plot(ax,[0 ztheta(3)], 'color', 'k');
        plot(ax,ztheta(3),'o','color','k','MarkerFaceColor','w', 'MarkerSize',14);
         plot(ax,[0 ztheta(4)], 'color', 'r');
        plot(ax,ztheta(4),'o','color','r','MarkerFaceColor','r', 'MarkerSize',14);
         plot(ax,[0 ztheta(5)], 'color', 'c');
        plot(ax,ztheta(5),'o','color','k','MarkerFaceColor','c', 'MarkerSize',14);
         plot(ax,[0 ztheta(6)], 'color', 'm');
        plot(ax,ztheta(6),'o','color','k','MarkerFaceColor','m', 'MarkerSize',14);
         plot(ax,[0 ztheta(7)], 'color', 'k');
        plot(ax,ztheta(7),'o','color','k','MarkerFaceColor','k', 'MarkerSize',14);
         plot(ax,[0 ztheta(8)], 'color', 'k');
        plot(ax,ztheta(8),'o','color','k','MarkerFaceColor','y', 'MarkerSize',14);
        
    else
        plot(ax,ztheta,'o','color','k');
    end
    % plot the centroid (gray paddle)
    plot(ax,[0 centroid], 'color', 'k');
    plot(ax,centroid,'o','color','k', 'Marker','o', 'MarkerFaceColor','k', 'MarkerSize',18);
    
    % axis limits etc
    axis(ax,'equal');
    xlim(ax,[-1.1 1.1]);
    ylim(ax,[-1.1 1.1]);
    % Added to be able to see the evolution of the phases
%     drawnow;
    %pause(0.01);
    %
end

% Auxiliary function for plotting the Kuramoto order parameter (R).
function KuramotoR(ax,t,sol,Kij,omega,Alpha)
    
% Project the phases into the complex plane.
    N_communities = 8;
    N_oscillators = 32;

    % Get the phases of the oscillators at time t
    theta = sol.y;
    ztheta = exp(1i.*theta);

    % Get the centroid for each community
    for c=1:N_communities
        if c==1
            ctheta(c,:)=mean(ztheta(1:32,:));
        else
            ctheta(c,:)=mean(ztheta((c-1)*N_oscillators:c*N_oscillators,:));

        end
    end
   

    % compute the running phase centroid
    centroid = ctheta;

    % plot the amplitide of the centroid versus time.
    axis(ax,'normal');
    plot(ax,sol.x,abs(centroid(1,:)),'color','b','linewidth',1);
    plot(ax,sol.x,abs(centroid(2,:)),'color','g','linewidth',1);
    plot(ax,sol.x,abs(centroid(3,:)),'color','k','linewidth',1);
    plot(ax,sol.x,abs(centroid(4,:)),'color','r','linewidth',1);
    plot(ax,sol.x,abs(centroid(5,:)),'color','c','linewidth',1);
    plot(ax,sol.x,abs(centroid(6,:)),'color','m','linewidth',1);
    plot(ax,sol.x,abs(centroid(7,:)),'color','k','linewidth',1);
    plot(ax,sol.x,abs(centroid(8,:)),'color','y','linewidth',1);
    
    % axis limits etc
    t0 = min(sol.x([1 end]));
    t1 = max(sol.x([1 end]));
    xlim(ax,[t0 t1]);
    ylim(ax,[-0.1 1.1]);
    xlabel(ax,'time');
    ylabel(ax,'R = abs(centroid)');
    title(ax,'Kuramoto Order Parameters for 8 oscillators');

    %hold on
end
