function Plot_phases_Shanahan_extra(start,stop,steps,record)
% function to plot the phase evolutions from the bdt simulation
%
% Fran Hancock
% April 2023
% fran.hancock@kcl.ac.uk
%
% Usage  Plot_phases_Shanahan_extra(1,2000,2,1,1)
%
%
PROJECT = '/Users/HDF/Hancock Dean Dropbox/Doug Dean/Fran/Academics/PostDoc/Projects/META_MODELLING';
addpath(genpath('/Users/HDF/Hancock Dean Dropbox/Doug Dean/Fran/Academics/PhD/Matlab_World/bdtoolbox-2022b'))

load('KuramotoShan.mat');

theta=sol.y;
N_communities = 8;
N_oscillators = 32;

t0 = start;
t1 = stop;

nrows = 6+3;

model = 2; %1 Hansel 2 Shanahan;
cmap = get_cmap(model);

if(record)
    videoModes = VideoWriter('videos/ShanahanChimera.mp4','MPEG-4'); %create the video object
    videoModes.FrameRate = 15;
    videoModes.Quality = 25;
    open(videoModes); %open the file for writing
end

hf=figure;
hf=colordef(hf,'white');
hf.Color='w';

set(gcf, 'units','normalized','outerposition',[0 0 0.25 1])

colormap jet


zztheta=exp(1i.*theta);

% get the centroids for each community

for c=1:N_communities
    if c==1
        ctheta(c,:)=mean(zztheta(1:32,:),1);
        ptheta(c,:) = angle(ctheta(c,:));

    else
        ctheta(c,:)=mean(zztheta((c-1)*N_oscillators:c*N_oscillators,:),1);
        ptheta(c,:) = angle(ctheta(c,:));
    end
end

for c=1:N_communities
    metakop_all(c,:) = var(abs(ctheta(c,:))); % KOP
end

ztheta= ctheta;
% compute the running phase centroid
centroid = ctheta;

for t=t0:steps:t1
    
    pos = 1; 

    subplot(nrows,1,[pos pos+2]);
    pos = pos + 2 + 2;

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
        plot([0 ztheta(1,t)], 'color', cmap{1});
        plot(ztheta(1,t),'o','color',cmap{1},'MarkerFaceColor',cmap{1}, 'MarkerSize',12);
        plot([0 ztheta(2,t)], 'color', cmap{2});
        plot(ztheta(2,t),'o','color',cmap{2},'MarkerFaceColor',cmap{2}, 'MarkerSize',12);
        plot([0 ztheta(3,t)], 'color', cmap{3});
        plot(ztheta(3,t),'o','color',cmap{3},'MarkerFaceColor',cmap{3}, 'MarkerSize',12);
        plot([0 ztheta(4,t)], 'color', cmap{4});
        plot(ztheta(4,t),'o','color',cmap{4},'MarkerFaceColor',cmap{4}, 'MarkerSize',12);
        plot([0 ztheta(5,t)], 'color', cmap{5});
        plot(ztheta(5,t),'o','color',cmap{5},'MarkerFaceColor',cmap{5}, 'MarkerSize',12);
        plot([0 ztheta(6,t)], 'color', cmap{6});
        plot(ztheta(6,t),'o','color',cmap{6},'MarkerFaceColor',cmap{6}, 'MarkerSize',12);
        plot([0 ztheta(7,t)], 'color', cmap{7});
        plot(ztheta(7,t),'o','color',cmap{7},'MarkerFaceColor',cmap{7}, 'MarkerSize',12);
        plot([0 ztheta(8,t)], 'color', cmap{8});
        plot(ztheta(8,t),'o','color',cmap{8},'MarkerFaceColor',cmap{8}, 'MarkerSize',12);
    else
        plot(ztheta,'o','color','k');
    end    
  
    % axis limits etc
    axis('equal');
    axis off
    xlim([-1.1 1.1]);
    ylim([-1.1 1.1]);

    txt=sprintf('Shanahan 2010 - Average phase for 8 communities of 32 oscillators \n \n t = %.2f secs   Brain Dynamics Toolbox v2022b',sol.x(t)); % 0.05 is th eintegration time steps)
    sgtitle(txt,'FontSize',16, 'FontWeight','b')

    
    ax = subplot(nrows,1, [pos pos+1]); % try to plot the evolution of KOP 
    pos = pos + 2 + 1;
    
    % plot the amplitide of the centroid versus time.
    axis('normal');
   
    cla(ax,"reset")

    hold on
    plot(ax,sol.x(t0:t1),abs(centroid(1,t0:t1)),'color',cmap{1},'linewidth',2);
    plot(ax,sol.x(t0:t1),abs(centroid(2,t0:t1)),'color',cmap{2},'linewidth',2);
    plot(ax,sol.x(t0:t1),abs(centroid(3,t0:t1)),'color',cmap{3},'linewidth',2);
    plot(ax,sol.x(t0:t1),abs(centroid(4,t0:t1)),'color',cmap{4},'linewidth',2);
    plot(ax,sol.x(t0:t1),abs(centroid(5,t0:t1)),'color',cmap{5},'linewidth',2);
    plot(ax,sol.x(t0:t1),abs(centroid(6,t0:t1)),'color',cmap{6},'linewidth',2);
    plot(ax,sol.x(t0:t1),abs(centroid(7,t0:t1)),'color',cmap{7},'linewidth',2);
    plot(ax,sol.x(t0:t1),abs(centroid(8,t0:t1)),'color',cmap{8},'linewidth',2);

    yyaxis(ax,"right")
    ymax = max(metakop_all,[],'all');
    ylim([0 ymax])
    ax1 = gca;
    ax1.YColor = 'k';
    ylabel('META_{KOP}')

    for mode=1:N_communities
        yline(metakop_all(mode,:),'--', 'LineWidth',2,'color',cmap{mode}) 
    end

    yyaxis(ax,"left")
    
    kt0 = min(sol.x([t0 t1])); % time in seconds
    kt1 = max(sol.x([t0 t1])); % time in seconds
    xlim(ax,[kt0 kt1]); % limits set wrt iterations
    xl = xlim;

    step1  = round((xl(2)-xl(1))/10); % used for ajusting the xaxis

    xticks([kt0:step1:kt1])
    xt = xticks;
    xlabs=string(round(xt));
    xticklabels(ax,xlabs)

    xline(ax,sol.x(t),'-', 'LineWidth',1,'color','r')

    xlabel(ax,'time');
    ylabel(ax,'KOP');
    title('Instantaneous Synchrony','FontSize',14,'FontWeight','bold');


    ax3=subplot(nrows,1, [pos pos+1]);
    cla
    plot_freq_sim('KuramotoShan.mat',start,stop,t,ax3,cmap,ptheta,8,32)
    pos = pos + 2 + 1;
        
 
    pause(0.001)


    if(record)
        frame = getframe(gcf);
        writeVideo(videoModes,frame); %write the image to file
    end
   
end



