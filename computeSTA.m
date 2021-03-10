function [sta] = computeSTA(spiketimes, stim, t_steps, DRAWNOW, GIF)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    f = logspace(log10(500), log10(22400), 34);
    stim_dim = size(stim, 1); %dimensionality of stimulus (i.e. 34 frequencies)
    
    n_trials = length(spiketimes);
    sta = zeros(stim_dim, t_steps);
    t_ax = - t_steps*6.25:6.25:0
    
    all_spike_times = [spiketimes{:}];
    n_spikes = length(all_spike_times);
    
    fig = figure;
    
    
    for i = 1:n_spikes
        ts = all_spike_times(i);
        col = ceil(ts/6.25);
        sta = sta + stim(:,col-t_steps+1:col);
        
        if DRAWNOW == true
            imagesc(t_ax, f, sta/i) %need to divide by i to normalise average
            cb = colorbar;
            %cb.Title.String('Average amplitude /dB')
            set(get(cb,'Title'),'String','Average amplitude /dB')
            axis xy
            title('Spike triggered average')
            xlabel('Time preceding spike /ms')
            ylabel('Frequency /Hz')
            set(gca, 'Fontsize', 14)

            drawnow;
        end
        
        %generates gif
        if GIF == true
            frame = getframe(fig); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            % Write to the GIF File 
            if i == 1 
              imwrite(imind,cm,'figs/STA.gif','gif', 'Loopcount',inf); 
            else 
              imwrite(imind,cm,'figs/STA.gif','gif','WriteMode','append'); 
            end 
        end
    end
    
    if DRAWNOW == false
        sta = sta/n_spikes;
        imagesc(t_ax, f, sta) %need to divide by i to normalise average
        cb = colorbar;
        %cb.Title.String('Average amplitude /dB')
        set(get(cb,'Title'),'String','Average amplitude /dB')
        axis xy
        title('Spike triggered average')
        xlabel('Time preceding spike /ms')
        ylabel('Frequency /Hz')
        set(gca, 'Fontsize', 14)
    end
end

