% Haken-Kelso-Bunz Model - Network of Phase Oscillators with 2nd order mode
%   Constructs a Kuramoto network with n nodes.
%       theta_i' = omega_i + SUM_j Kij*sin(theta_i-theta_j) - SUM_j
%       Kij*sin(2(theta_i-theta_j)
%   where 
%       theta is an (nx1) vector of oscillator phases (in radians),
%       omega is an (nx1) vector of natural frequencies (cycles/sec)
%       Kij is an (nxn) matrix of connection weights,
%
% Example:
%   n = 8;                    % number of oscillators
%   Kij = ones(n);             % coupling matrix
%   sys = HKB(Kij);    % construct the system struct
%   gui = bdGUI(sys);          % open the Brain Dynamics GUI
%
% Usage
% For a=b=0.105 Zhang et al (2019)
% n=8;Kij=ones(n);Kij=Kij-diag(diag(Kij)); Kij=0.105*Kij;sys=HKB(Kij);gui=bdGUI(sys)
% For a=b=0.15 Zhang thesis
% n=8;Kij=ones(n);Kij=Kij-diag(diag(Kij)); Kij=0.15*Kij;sys=HKB(Kij);gui=bdGUI(sys);
%
% Author: Fran Hancock (2023)
% Adapted from:
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
function sys = HKB(Kij)
    % Determine the number of nodes from the size of the coupling matrix
    addpath(genpath('/Users/HDF/Hancock Dean Dropbox/Doug Dean/Fran/Academics/PhD/Matlab_World/bdtoolbox-2022b'))

    n = size(Kij,1);
    dw = 0.075;
    w=1.5;
    
    % Handle to our ODE function
    sys.odefun = @odefun;
    
    % ODE parameters
    sys.pardef = [ struct('name','Kij',   'value',Kij);
                   struct('name','omega', 'value',[w-3*dw; w-2*dw; w-dw; w; w+dw; w+2*dw; w+3*dw; w]);
                   ];
               
    % ODE state variables
    sys.vardef = struct('name','theta', 'value',2*pi*rand(n,1), 'lim',[0 2*pi]);
    
    % Time span
    sys.tspan = [0 1000];
    sys.odeoption.MaxStep = 0.04; 
    sys.tstep = 0.04;

    % Relevant ODE solvers
    sys.odesolver = {@ode45,@ode23,@ode113,@odeEul};
    
    % ODE solver options
    sys.odeoption.RelTol = 1e-6;
    %sys.odeoption.MaxStep = 0.1;                
                    
    % Latex (Equations) panel
    sys.panels.bdLatexPanel.title = 'Equations'; 
    sys.panels.bdLatexPanel.latex = {
        '$\textbf{Generalized HKB}$'
        ''
        'A generalised network of HKB Oscillators with 2nd order coupling'
        '{ }{ }{ } $\dot \theta_i = \omega_i - a\sum_j K_{ij} \sin(\theta_i - \theta_j) - b\sum_j K_{ij} \sin(2*(\theta_i - \theta_j))$'
        'where'
        '{ }{ }{ } $\theta_i$ is the phase of the $i^{th}$ oscillator (radians),'
        '{ }{ }{ } $\omega_i$ is a distribution of oscillator frequencies (cycles/sec),'
        '{ }{ }{ } $i,j=1 \dots n,$'
        ['{ }{ }{ } $n{=}',num2str(n),'.$']
        ''
        'References'
        'Kuramoto (1984) Chemical oscillations, waves and turbulence.'
        'Strogatz (2000) From Kuramoto to Crawford.'
        'Zjang et al (2019) Connecting empirical phenomena and theoretical models of biological coordination across scales.'
        'Chapter 6.2 of the Handbook for the Brain Dynamics Toolbox (Version 2018b).'
        };
    
    % Time Portrait panel
    sys.panels.bdTimePortrait.title = 'Time Portrait';
    sys.panels.bdTimePortrait.modulo = 'on';
 
    % Phase Portrait panel
    sys.panels.bdPhasePortrait.title = 'Phase Portrait';
    sys.panels.bdPhasePortrait.modulo = 'on';

    % Auxiliary panel
    sys.panels.bdAuxiliary.title = 'Auxiliary';
    sys.panels.bdAuxiliary.auxfun = {@centroid1,@centroid2,@KuramotoR,@Frequency};


    % Solver panel
    sys.panels.bdSolverPanel.title = 'Solver';                
end

% Kuramoto ODE function where
% theta is a (1xn) vector of oscillator phases,
% Kij is either a scalar or an (nxn) matrix of connection weights,
% omega is either a scalar or (1xn) vector of oscillator frequencies.

function dtheta = odefun(t,theta,Kij,omega)
    n = numel(theta);
    theta_i = theta * ones(1,n);                        % (nxn) matrix with same theta values in each row
    theta_j = ones(n,1) * theta';                       % (nxn) matrix with same theta values in each col
    theta_ij = theta_i - theta_j;                       % (nxn) matrix of all possible (theta_i - theta_j) combinations
    dtheta = omega - sum(Kij.*sin(theta_ij),1)' - sum(Kij.*sin(2.*(theta_ij)),1)' ;   % Kuramoto Equation in vector form.

    %  dtheta = omega - sum(a.*Kij.*sin(theta_ij),1)' - sum(b.*Kij.*sin(2*(theta_ij)),1)';   % Kuramoto Equation in vector form.
end

% Auxiliary function that plots the centroid of the oscillators
function centroid1(ax,t,sol,Kij,omega)

    % Get the phases of the oscillators at time t
    theta = bdEval(sol,t);
    
    % Project the phases into the complex plane.
    ztheta = exp(1i.*theta);
    
    % Plot the centroid.
    centroidplot(ax,ztheta);
    %text(ax,-1,-1,num2str(t,'time = %g'));
    title(ax,'centroid of oscillators'); 
end

% Auxiliary function that plots the centroid of the oscillators
% in a rotating frame where the first oscillator is pinned at
% zero phase.
function centroid2(ax,t,sol,Kij,omega)
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
    title(ax,'centroid (rotating frame)'); 
end

% Utility function for plotting the centroid
function centroidplot(ax,ztheta)
    
    cmap = get_cmap;
% compute the phase centroid
    centroid = mean(ztheta);
    
    % plot the unit circle
    plot(ax,exp(1i.*linspace(-pi,pi,100)), 'color',[0.75 0.75 0.75]);
    
    % plot the oscillator phases on the unit circle
    plot(ax,[0 ztheta(1)], 'color', cmap{1});
    plot(ax,ztheta(1),'o','color',cmap{1},'MarkerFaceColor',cmap{1}, 'MarkerSize',12);
    
    plot(ax,[0 ztheta(2)], 'color', cmap{2});
    plot(ax,ztheta(2),'o','color',cmap{2},'MarkerFaceColor',cmap{2}, 'MarkerSize',12);
    
    plot(ax,[0 ztheta(3)], 'color', cmap{3});
    plot(ax,ztheta(3),'o','color',cmap{3},'MarkerFaceColor',cmap{3}, 'MarkerSize',12);
    
    plot(ax,[0 ztheta(4)], 'color', cmap{4});
    plot(ax,ztheta(4),'o','color',cmap{4},'MarkerFaceColor',cmap{4}, 'MarkerSize',12);
    
    plot(ax,[0 ztheta(5)], 'color', cmap{5});
    plot(ax,ztheta(5),'o','color',cmap{5},'MarkerFaceColor',cmap{5}, 'MarkerSize',12);
    
    plot(ax,[0 ztheta(6)], 'color', cmap{6});
    plot(ax,ztheta(6),'o','color',cmap{6},'MarkerFaceColor',cmap{6}, 'MarkerSize',12);
    
    plot(ax,[0 ztheta(7)], 'color', cmap{7});
    plot(ax,ztheta(7),'o','color',cmap{7},'MarkerFaceColor',cmap{7}, 'MarkerSize',12);
    
    plot(ax,[0 ztheta(8)], 'color', cmap{8});
    plot(ax,ztheta(8),'o','color',cmap{8},'MarkerFaceColor',cmap{8}, 'MarkerSize',12);

    % plot the centroid (yellow paddle)
%     plot(ax,[0 centroid], 'color', 'y');
%     plot(ax,centroid,'o','color','y', 'Marker','o', 'MarkerFaceColor','y', 'MarkerSize',10);
    
    % axis limits etc
    ylabel(ax,'');
    axis(ax,'equal');
    xlim(ax,[-1.1 1.1]);
    ylim(ax,[-1.1 1.1]);
    % Added to be able to see the evolution of the phases
  

    pause(0.01);
    %
end

% Auxiliary function for plotting the Kuramoto order parameter (R).
function KuramotoR(ax,t,sol,Kij,omega)
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
    title(ax,'Kuramoto Order Parameter (R1 and R2)');
    
    hold on

    ztheta = exp(1i.*2*sol.y);

    % compute the running phase centroid
    centroid = mean(ztheta);

    % plot the amplitide of the centroid versus time.
    axis(ax,'normal');
    plot(ax,sol.x,abs(centroid),'color','b','linewidth',1);
end

function Frequency(ax,t,sol,Kij,omega)
    % Project the phases into the complex plane.
    
    cmap = get_cmap;
    theta = sol.y';
    time = sol.x';

    freq = (diff(theta)./diff(time));
    freq = freq';
    time = time(2:end)';
   
    % plot the amplitide of the centroid versus time.
    axis(ax,'normal');
    plot(ax,time,freq(1,:),'color',cmap{1},'linewidth',1);
    plot(ax,time,freq(2,:),'color',cmap{2},'linewidth',1);
    plot(ax,time,freq(3,:),'color',cmap{3},'linewidth',2);
    plot(ax,time,freq(4,:),'color',cmap{4},'linewidth',1);
    plot(ax,time,freq(5,:),'color',cmap{5},'linewidth',1);
    plot(ax,time,freq(6,:),'color',cmap{6},'linewidth',1);
    plot(ax,time,freq(7,:),'color',cmap{7},'linewidth',1);
    plot(ax,time,freq(8,:),'color',cmap{8},'linewidth',1);

    % axis limits etc
    t0 = time(1);
    t1 = time(end);
    xlim(ax,[t0 t1]);
    
    ymin = min(freq,[],'all');
    if ymin == 0
        ymin = -0.1;
    end
    
    ymax = max(freq,[],'all');

    if ymax>0
        ylim(ax,[ymin ymax]);
    else
        ylim(ax,[ymin 1]);
    end
    xlabel(ax,'time seconds');
    ylabel(ax,'Frequency Hz');
    title(ax,'Instantaneous frequency');
end

% Auxiliary function for plotting the Kuramoto order parameter (R2).
function KuramotoCombined(ax,t,sol,Kij,omega)
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
    title(ax,'Kuramoto Order Parameter (R)');
end

function [cmap] = get_cmap
    blue = [0 0 255]./256;
    green = [0.4660 0.6740 0.1880];
    white = [200 200 200]./256; %grey
    red = [255 0 0]./256;
    cyan = [0 255 255]./256;
    magenta = [0.4940 0.1840 0.5560];
    black = [0 0 0];
    yellow = [255 255 0]./256;
    dk_blue = [0 0.4470 0.7410];
    orange = [0.9290 0.6940 0.1250];
    cherry = [0.6350 0.0780 0.1840];

    cmap= {blue, cyan, green, dk_blue, orange,red, cherry, magenta};
    
end