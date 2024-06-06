% Hansel-Mato-Meunier model - Network of Kuramoto Phase Oscillators with 2nd order mode
%   Constructs a Kuramoto network with n nodes.
%       theta_i' = omega_i + SUM_j Kij*sin(theta_i-theta_j+ Alpha) + SUM_j Kij*sin(2(theta_i-theta_j)) + Eta
%   where 
%       theta is an (nx1) vector of oscillator phases (in radians),
%       omega is an (nx1) vector of natural frequencies (cycles/sec)
%       Kij is an (nxn) matrix of connection weights,
%       Alpha is a fixed phase lag
%       Eta is a random noise term
%
% Example:
%   n = 8;                    % number of oscillators
%   Kij = ones(n);             % coupling matrix
%   sys = Kuramoto2order(Kij);    % construct the system struct
%   gui = bdGUI(sys);          % open the Brain Dynamics GUI
%   n=5; kij=rand(n,n); Kij=kij-diag(diag(kij)) + diag(1); sys=Kuramoto2order(Kij); gui=bdGUI(sys);
% 
% Author: Fran Hancock, 2023 
% Adapted from:
%
% Authors
%   Stewart Heitmann (2016a,2017a,2018a,2018b,2020a)

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
function sys = HMM(Kij)
    addpath(genpath('/Users/HDF/Hancock Dean Dropbox/Doug Dean/Fran/Academics/PhD/Matlab_World/bdtoolbox-2022b'))
    % Determine the number of nodes from the size of the coupling matrix
    n = size(Kij,1);
    
    % Handle to our ODE function
    sys.odefun = @odefun;
    
    % ODE parameters
    sys.pardef = [ struct('name','Kij',   'value',Kij);
               %    struct('name','k',     'value',1); 
                   struct('name','omega', 'value',1);%[1.1600; 3.6698; 0.6771; 5.6945; 5.527]); %randn(n,1));
                   struct('name','R', 'value',0.25); % as in Hansel et al 1993
               %    struct('name','Beta','value',0); % as in Hansel 1993
                   struct('name','Alpha', 'value',1.435); % as in Ashwin 2007 - 3 cluster states
                   struct('name','Eta', 'value',sqrt(0.0022).*randn(n,1)); % as in Ashwin 2007 - 3 cluster states
                   ];
               
    % ODE state variables
    sys.vardef = struct('name','theta', 'value',[1.1600; 3.6698; 0.6771; 5.6945; 5.527]);
   %sys.vardef = struct('name','theta', 'value',2*pi*rand(n,1), 'lim',[0 2*pi]);
    
%     sys.tspan = [0 100];
%     sys.tstep = 0.03;
    sys.tspan = [0 1000];
    sys.tstep = 0.05;
    
    % Relevant ODE solvers
    sys.odesolver = {@ode45,@ode23,@ode113,@odeEul};
    
    % ODE solver options
    sys.odeoption.RelTol = 1e-6;
    sys.odeoption.MaxStep = 0.03;
                   
                    
    % Latex (Equations) panel
    sys.panels.bdLatexPanel.title = 'Equations'; 
    sys.panels.bdLatexPanel.latex = {
        '$\textbf{Kuramoto2order}$'
        ''
        'A generalised network of Kuramoto Oscillators with 2nd order coupling'
        '{ }{ }{ } $\dot \theta_i = \omega_i - \frac{k}{n} \sum_j K_{ij} \sin(\theta_i - \theta_j + Alpha) + \frac{k}{n} \sum_j K_{ij} \sin(2*(\theta_i - \theta_j)) + eta $'
        'where'
        '{ }{ }{ } $\theta_i$ is the phase of the $i^{th}$ oscillator (radians),'
        '{ }{ }{ } $\omega_i$ is its natural oscillation frequency (cycles/sec),'
        '{ }{ }{ } $K$ is the network connectivity matrix ($n$ x $n$),'
        '{ }{ }{ } $Aplha$ is a constant phase lag,'
        '{ }{ }{ } $Eta$ is random noise,'
        '{ }{ }{ } $i,j=1 \dots n,$'
        ['{ }{ }{ } $n{=}',num2str(n),'.$']
        ''
        'References'
        'Kuramoto (1984) Chemical oscillations, waves and turbulence.'
        'Hansel et al (1993) Clustering and slow switching in globally coupled phase oscillators'
        'Breakspear et al (2010) Generative models of cortical oscillations.'
        'Chapter 6.2 of the Handbook for the Brain Dynamics Toolbox (Version 2018b).'
        'Ashwin et al (2007) Dynamics on Networks of Cluster States for Globally Coupled Phase Oscillators.'
        };
    
    % Time Portrait panel
    sys.panels.bdTimePortrait.title = 'Time Portrait';
    sys.panels.bdTimePortrait.modulo = 'on';
 
    % Phase Portrait panel
    sys.panels.bdPhasePortrait.title = 'Phase Portrait';
    sys.panels.bdPhasePortrait.modulo = 'on';

    % Auxiliary panel
    sys.panels.bdAuxiliary.title = 'Auxiliary';
    sys.panels.bdAuxiliary.auxfun = {@centroid1,@centroid2,@KuramotoR,@KuramotoR2,@KuramotoCombined};


    % Solver panel
    sys.panels.bdSolverPanel.title = 'Solver';                
end

% Kuramoto ODE function where
% theta is a (1xn) vector of oscillator phases,
% Kij is either a scalar or an (nxn) matrix of connection weights,
% k is a scalar,
% omega is either a scalar or (1xn) vector of oscillator frequencies.
function dtheta = odefun(t,theta,Kij,omega,R, Alpha, Eta)
    n = numel(theta);
    theta_i = theta * ones(1,n);                        % (nxn) matrix with same theta values in each row
    theta_j = ones(n,1) * theta';                       % (nxn) matrix with same theta values in each col
    theta_ij = theta_i - theta_j;                       % (nxn) matrix of all possible (theta_i - theta_j) combinations
    dtheta = omega - 1/n.*sum(Kij.*sin(theta_ij+Alpha),1)' + 1/n.*R*sum(Kij.*sin(2*(theta_ij)),1)' + Eta;   % Kuramoto Equation in vector form.

end

% Auxiliary function that plots the centroid of the oscillators
function centroid1(ax,t,sol,Kij,omega,R,Alpha,Eta)

    % Get the phases of the oscillators at time t
    theta = bdEval(sol,t);
    
    % Project the phases into the complex plane.
    ztheta = exp(1i.*theta);
    
    % Plot the centroid.
    centroidplot(ax,ztheta);
    %text(ax,-1,-1,num2str(t,'time = %g'));
    title(ax,'centroid of oscillators'); 
end

% % Auxiliary function that plots the centroid of the oscillators
% % in a rotating frame where the first oscillator is pinned at
% % zero phase.
function centroid2(ax,t,sol,Kij,omega,R,Alpha,Eta)
    % Get the phases of the oscillators at time t
    theta = bdEval(sol,t);
    
    % Project the phases into the complex plane and
    % rotate the frame to pin the first oscillator at zero phase.
    ztheta = exp(1i.*theta) .* exp(-1i*theta(1));
    [wrow, ~]=find(angle(ztheta)==0);
    if wrow > 0
        ztheta(wrow,:)=[1 + 0.01i]; % hack to get the line to stay on the unit circle
    end
        % Plot the centroid.
    centroidplot(ax,ztheta);
    text(ax,-1,-1,num2str(t,'time = %g'));
    title(ax,'centroid (Relative Phase)'); 
end

function centroid3(ax,t,sol,Kij,omega,R,Alpha,Eta)
    % Plot relative phase
    % Get the phases of the oscillators at time t
    theta = bdEval(sol,t);
    
    % Project the phases into the complex plane and
    rtheta=theta;
    p0=rtheta(1,1);
    rtheta=rtheta - rtheta(1,:);

    ztheta = exp(1i.*rtheta);
    ztheta(1,1)=[1 + 0.01i]; % hack to get the line to stay on the unit circle
    % now calculate the relative phases - code taken from Hilbert function

%     %p0=ztheta(1,1);
%     rtheta=unwrap(ztheta - ztheta(1,:)); % +p0
    
    % Plot the centroid.
    centroidplot(ax,ztheta);
    text(ax,-1,-1,num2str(t,'time = %g'));
    title(ax,'centroid (relative phase)'); 
end

% Utility function for plotting the centroid
function centroidplot(ax,ztheta)
    % compute the phase centroid
    centroid = mean(ztheta);
    
    % plot the unit circle
    plot(ax,exp(1i.*linspace(-pi,pi,100)), 'color',[0.75 0.75 0.75]);
    
    % plot the oscillator phases on the unit circle
    if size(ztheta,1)==5
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
    else
        plot(ax,ztheta,'o','color','k');
    end
    % plot the centroid (yellow paddle)
    plot(ax,[0 centroid], 'color', 'k');
    plot(ax,centroid,'o','color','k', 'Marker','o', 'MarkerFaceColor','k', 'MarkerSize',14);
    
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
function KuramotoR(ax,t,sol,Kij,omega,R,Alpha,Eta)
    % Project the phases into the complex plane.
    ztheta = exp(1i.*sol.y);

    % compute the running phase centroid
    centroid = mean(ztheta);

    % plot the amplitide of the centroid versus time.
    axis(ax,'normal');
    plot(ax,sol.x,abs(centroid),'color','r','linewidth',1);
    
    % axis limits etc
    t0 = min(sol.x([1 end]));
    t1 = max(sol.x([1 end]));
    xlim(ax,[t0 t1]);
    ylim(ax,[-0.1 1.1]);
    xlabel(ax,'time');
    ylabel(ax,'R = abs(centroid)');
    title(ax,'Kuramoto Order Parameters (R1 R2)');

    %hold on
    
    ztheta = exp(1i.*2*sol.y);

    % compute the running phase centroid
    centroid = mean(ztheta);

    % plot the amplitide of the centroid versus time.
    %axis(ax,'normal');
    plot(ax,sol.x,abs(centroid),'color','b','linewidth',1);
    

end

% Auxiliary function for plotting the Kuramoto order parameter (R2).
function KuramotoR2(ax,t,sol,Kij,omega,R,Alpha,Eta)
    % Project the phases into the complex plane.
    ztheta = exp(1i.*2*sol.y);

    % compute the running phase centroid
    centroid = mean(ztheta);

    % plot the amplitide of the centroid versus time.
    axis(ax,'normal');
    plot(ax,sol.x,abs(centroid),'color','k','linewidth',1);
    
    % axis limits etc
    t0 = min(sol.x([1 end]));
    t1 = max(sol.x([1 end]));
    xlim(ax,[t0 t1]);
    ylim(ax,[-0.1 1.1]);
    xlabel(ax,'time');
    ylabel(ax,'R = abs(centroid)');
    title(ax,'Kuramoto Order Parameter (R2)');
end

% Auxiliary function for plotting the Kuramoto order parameter (R2).
function KuramotoCombined(ax,t,sol,Kij,omega,R,Alpha,Eta)
    % Project the phases into the complex plane.
    
    z1theta = exp(1i.*sol.y);
    z2theta = exp(1i.*2*sol.y);
    ztheta=(z1theta+z2theta)/2;
    % compute the running phase centroid
    centroid = mean(ztheta);

    % plot the amplitide of the centroid versus time.
    axis(ax,'normal');
    plot(ax,sol.x,abs(centroid),'color','k','linewidth',1);
    
    % axis limits etc
    t0 = min(sol.x([1 end]));
    t1 = max(sol.x([1 end]));
    xlim(ax,[t0 t1]);
    ylim(ax,[-0.1 1.1]);
    xlabel(ax,'time');
    ylabel(ax,'R = abs(centroid)');
    title(ax,'Kuramoto Order Parameter ((R1 + R2) /2)');
end