function plot_freq_sim(sim_results,start,stop,tnow,ax,cmap,theta,N_communities, N_oscillators)

load(sim_results)
t0 = start;
t1 = stop;

theta = sol.y(:,start:stop)';

% theta = sol.y(:,t0:t1)';
time = sol.x(1,t0:t1)';

freq = abs((diff(theta)./diff(time)));
freq = freq';

if size(theta,2) == 256 % Shanahan model

    % Calculate the average frequency for each community
    for c = 1:N_communities
        if c == 1
            cfreq(c,:) = mean(freq(1:N_oscillators,:));
        else
            cfreq(c,:) = mean(freq((c-1)*N_oscillators:c*N_oscillators,:));
        end
    end
    freq = cfreq;
end

time = time(2:end)';



% plot the amplitide of the centroid versus time.
axis(ax,'normal');

hold on
if size(freq,1) == 5 % Hansel
    plot(ax,time,freq(1,:),'color',cmap{1},'linewidth',2);
    plot(ax,time,freq(2,:),'color',cmap{2},'linewidth',2);
    plot(ax,time,freq(3,:),'color',cmap{3},'linewidth',2);
    plot(ax,time,freq(4,:),'color',cmap{4},'linewidth',2);
    plot(ax,time,freq(5,:),'color',cmap{5},'linewidth',2);

elseif size(freq,1) == 8 
    plot(ax,time,freq(1,:),'color',cmap{1},'linewidth',2);
    plot(ax,time,freq(2,:),'color',cmap{2},'linewidth',2);
    plot(ax,time,freq(3,:),'color',cmap{3},'linewidth',2);
    plot(ax,time,freq(4,:),'color',cmap{4},'linewidth',2);
    plot(ax,time,freq(5,:),'color',cmap{5},'linewidth',2);
    plot(ax,time,freq(6,:),'color',cmap{6},'linewidth',2);
    plot(ax,time,freq(7,:),'color',cmap{7},'linewidth',2);
    plot(ax,time,freq(8,:),'color',cmap{8},'linewidth',2);

end

% axis limits etc

%axis limits etc
t0 = time(1);
t1 = time(end);
xlim(ax,[t0 t1]);

if size(theta,2) == 8
    ymax = 1.8; 
    ymin = 1.2;
else
    ymax = max(freq,[],'all');
    ymin = min(freq,[],'all');
    if ymin == 0
        ymin = -0.1;
    end
end
if ymax>0
    ylim(ax,[ymin ymax]);
else
    ylim(ax,[ymin 1]);
end

xline(sol.x(tnow),'-', 'LineWidth',1,'color','r')

xlabel(ax,'time in seconds', 'FontSize',16);
ylabel(ax,'Frequency Hz','FontSize',16);
title(ax,['Instantaneous frequency' newline],'FontSize',16);
