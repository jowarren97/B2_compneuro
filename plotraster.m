function [] = plotraster(spiketimes)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    if class(spiketimes) == "cell"
        n_trials = length(spiketimes);
        figure('Position', [100, 100, 1000, n_trials*25])
        hold on
        for i = 1:n_trials %loop thru trials
            indiv_trial = spiketimes{i}; %select data of trial
            sz = size(indiv_trial); %find length of trial
            plot(spiketimes{i}, i*ones(sz), 'ko', 'MarkerSize', 4)
        end
        hold off
        
    elseif class(spiketimes) == "double" %i.e. if single trial
        n_trials = 1;
        figure('Position', [100, 100, 1000, 200])
        sz = size(spiketimes);
        plot(spiketimes, ones(sz), 'ko', 'MarkerSize', 4)
    end   
    
    xlabel('time /ms')
    ylabel('Trial number')
    title('Raster plot of spike times')
    ylim([0.5 n_trials+0.5])
    set(gca,'Ytick',1:1:n_trials)
    set(gca,'YtickLabel',1:1:n_trials)
    set(gca, 'Fontsize', 14)
end

