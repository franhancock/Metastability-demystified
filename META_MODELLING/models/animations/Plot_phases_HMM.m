function Plot_phases_HMM(start,stop,steps,record)
% function to plot the phase evolutions from the bdt simulation
%
% Fran Hancock
% April 2023
% fran.hancock@kcl.ac.uk
%
% Usage Plot_phases_HMM(1,30001,100,1)
% Use lower time steps for shorter intervals eg (1,20001,50,1)
%

addpath(genpath('/Users/HDF/Hancock Dean Dropbox/Doug Dean/Fran/Academics/PhD/Matlab_World/bdtoolbox-2022b'))

t0 = start;
t1 = stop;
model = 1; %1 Hansel 2 Shanahan;
cmap = get_cmap(model);

%if(~reduced)    
%    nrows = 15;
%else
    nrows = 6;
%end
pos = 1; % allow differnt figures

load HMM.mat
theta = sol.y;

if(record)

    videoModes = VideoWriter('videos/HMM.mp4','MPEG-4'); %create the video object
    videoModes.FrameRate = 15;
    videoModes.Quality = 25;
    open(videoModes); %open the file for writing
end

hf=figure;
hf=colordef(hf,'white');
hf.Color='w';

set(gcf, 'units','normalized','outerposition',[0 0 0.25 0.5])

for t=t0:steps:t1

    pos = 1; 

    subplot(nrows,1,[pos pos+2])
    pos = pos + 2 + 2;

    cla
    hold on

    %ztheta = (exp(1i.*theta(:,t)) .* exp(-1i*theta(1,t))); % relative phase
    ztheta = exp(1i*theta(:,t));
    [wrow, ~]=find(angle(ztheta)==0);
    if wrow > 0
        ztheta(wrow,:)=[1 + 0.01i]; % hack to get the line to stay on the unit circle
    end
    
    cla
    hold on

    plot(exp(1i.*linspace(-pi,pi,100)), 'color','k'); % plot the outer circle
    
    th = 0:pi/50:2*pi; % Plot the inner circel
    x_circle = 0.5 * cos(th) ;
    y_circle = 0.5 * sin(th);
    plot(x_circle, y_circle, 'color','k','LineStyle','--');

    plot([0 0+1i], 'color','k','LineStyle','--')
    plot([0 0-1i], 'color','k','LineStyle','--')
    plot([0 1+0.01i], 'color','k','LineStyle','--')
    plot([0 -1+0.01i], 'color','k','LineStyle','--')
    plot([0 sqrt(3)/2+0.5i], 'color','k','LineStyle','--')
    plot([0 0.5+sqrt(3)/2*i], 'color','k','LineStyle','--')
    plot([0 -sqrt(3)/2+0.5i], 'color','k','LineStyle','--')
    plot([0 -0.5+sqrt(3)/2*i], 'color','k','LineStyle','--')
    plot([0 sqrt(3)/2-0.5i], 'color','k','LineStyle','--')
    plot([0 0.5-sqrt(3)/2*i], 'color','k','LineStyle','--')
    plot([0 -sqrt(3)/2-0.5i], 'color','k','LineStyle','--')
    plot([0 -0.5-sqrt(3)/2*i], 'color','k','LineStyle','--')
   

    % plot the oscillator phases on the unit circle
    if size(ztheta,1)==5
        plot([0 ztheta(1)], 'color', cmap{1},'LineWidth',1.5);
        plot(ztheta(1),'o','color',cmap{1},'MarkerFaceColor',cmap{1}, 'MarkerSize',14);
         plot([0 ztheta(2)], 'color', cmap{2},'LineWidth',1.5);
        plot(ztheta(2),'o','color',cmap{2},'MarkerFaceColor',cmap{2}, 'MarkerSize',14);
         plot([0 ztheta(3)], 'color', cmap{3},'LineWidth',1.5);
        plot(ztheta(3),'o','color',cmap{3},'MarkerFaceColor',cmap{3}, 'MarkerSize',14);
         plot([0 ztheta(4)], 'color', cmap{4},'LineWidth',1.5);
        plot(ztheta(4),'o','color',cmap{4},'MarkerFaceColor',cmap{4}, 'MarkerSize',14);
         plot([0 ztheta(5)], 'color', cmap{5},'LineWidth',1.5);
        plot(ztheta(5),'o','color',cmap{5},'MarkerFaceColor',cmap{5}, 'MarkerSize',14);
    else
        plot(ztheta,'o','color','k');
    end
    
  
    % axis limits etc
    axis('equal');
    axis off
    xlim([-1.1 1.1]);
    ylim([-1.1 1.1]);
    
    txt=sprintf('Relative phase for Hansel-Mato-Meunier Model (1993) \n \n t =  = %.2f secs   Brain Dynamics Toolbox v2022b \n \n',sol.x(t)); % 0.05 is th eintegration time steps)
    sgtitle(txt,'FontSize',16, 'FontWeight','b')
    
    % FREQUENCY
    ax=subplot(nrows,1, [pos pos+1]);
    pos = pos + 2 + 1;
    cla
    
%    plot_freq_sim(sim_results,start,stop,tnow,ax,cmap,theta,N_communities, N_oscillators)
    plot_freq_sim('HMM.mat',start,stop,t,ax,cmap,theta)

    
%     if(~reduced)
%         % EIGENVECTORS
%         ax = subplot(nrows,1,[pos pos+1]);
%         pos = pos + 2 + 1;
%         cla
%         
%         if t == t0
%             Plot_mode_eigs_sim('HanselFreq.mat',start,stop,0,0,0.6,1,1,ax,cmap,seppe,tpose)
%         else
%             %function Plot_mode_eigs_sim(sim_results,start,stop,N_communities,N_oscillators,range,eigV,tnow,ax)
%             Plot_mode_eigs_sim('HanselFreq.mat',start,stop,0,0,0.6,0,t,ax,cmap,seppe,tpose)
%         end
%     end

    
%     if(~reduced)
%         % Phase Difference Differential
%         ax=subplot(nrows,1, [pos pos + 1]);
%         pos = pos + 2 + 1;
%         cla
%  
%         Plot_PDDBimodal('HanselFreq.mat',start,stop,t,ax,cmap,1,2)
%      
%         ax = subplot(nrows,1, [pos pos + 1]);
%         cla
%         Plot_PDD('HanselFreq.mat',start,stop,t,ax,cmap,0,0)
%         
%          if(seppe && tpose)
%          %plot_eigenvalues('HanselFreq.mat', ax2, start, stop,t)
%         end
%     end

    if(record)
        frame = getframe(gcf);
        writeVideo(videoModes,frame); %write the image to file
    end
     if t == t0
         pause(0.05)
     else
         pause(0.0001)
     end 
    
end

