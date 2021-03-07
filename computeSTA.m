function [] = computeSTA(spiketimes, stim, t_steps)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    f = logspace(log10(500), log10(22400), 34);
    stim_dim = size(stim, 1); %dimensionality of stimulus (i.e. 34 frequencies)
    
    n_trials = length(spiketimes);
    sta = zeros(stim_dim, t_steps);
    t_ax = - t_steps*6.25:6.25:0
    
    all_spike_times = [spiketimes{:}];
    n_spikes = length(all_spike_times);
    
    figure    
    for i = 1:n_spikes
        ts = all_spike_times(i);
        col = ceil(ts/6.25);
        sta = sta + stim(:,col-15:col-1);
        imagesc(t_ax, f, sta/i) %need to divide by i to normalise average
        drawnow;
    end
    
    cb = colorbar;
    cb.Title.String('Average amplitude /dB')
    axis xy
    title('Spike triggered average')
    xlabel('Time preceding spike /ms')
    ylabel('Frequency /Hz')
end

