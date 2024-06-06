function Plot_phases_HKB(start,stop,steps,record)
% function to plot the phase evolutions from the bdt simulation
%
% Fran Hancock
% April 2023
% fran.hancock@kcl.ac.uk
%
% Usage Plot_phases_HKB(1,6000,10,1,1,1,0) for video
% Use lower time steps for shorter intervals eg (1,20001,50)
%

addpath(genpath('/Users/HDF/Hancock Dean Dropbox/Doug Dean/Fran/Academics/PhD/Matlab_World/bdtoolbox-2022b'))

t0 = start;
t1 = stop;
model = 2; %1 Hansel Shanahan 2 HKB;
cmap = get_cmap(model);

nrows = 6;
pos = 1; % allow differnt figures

load HKB
theta = sol.y;

if(record)
    videoModes = VideoWriter('videos/HKB.mp4','MPEG-4'); %create the video object
    videoModes.FrameRate = 4; %15;
    videoModes.Quality = 25;
    open(videoModes); %open the file for writing
end

hf=figure;
hf=colordef(hf,'white');
hf.Color='w';

set(gcf, 'units','normalized','outerposition',[0 0 0.5 0.5])

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
    if size(ztheta,1)==8
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
         plot([0 ztheta(6)], 'color', cmap{6},'LineWidth',1.5);
        plot(ztheta(6),'o','color',cmap{6},'MarkerFaceColor',cmap{6}, 'MarkerSize',14);
         plot([0 ztheta(7)], 'color', cmap{7},'LineWidth',1.5);
        plot(ztheta(7),'o','color',cmap{7},'MarkerFaceColor',cmap{7}, 'MarkerSize',14);
         plot([0 ztheta(8)], 'color', cmap{8},'LineWidth',1.5);
        plot(ztheta(8),'o','color',cmap{8},'MarkerFaceColor',cmap{8}, 'MarkerSize',14);
    else
        plot(ztheta,'o','color','k');
    end
    
  
    % axis limits etc
    axis('equal');
    axis off
    xlim([-1.1 1.1]);
    ylim([-1.1 1.1]);
    
    txt=sprintf('Relative phase for Generalised HKB \n \n t = %.2f secs   Brain Dynamics Toolbox v2022b \n \n',sol.x(t)); % 0.05 is th eintegration time steps)
    sgtitle(txt,'FontSize',16, 'FontWeight','b')
    
    % FREQUENCY
    ax=subplot(nrows,1, [pos pos+1]);
    pos = pos + 2 + 1;
    cla
    plot_freq_sim('HKB.mat',start,stop,t,ax,cmap,theta)


    if(record)
        frame = getframe(gcf);
        writeVideo(videoModes,frame); %write the image to file
    end

    pause(0.05)
    
end

