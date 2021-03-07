function [] = plotraster(spiketimes)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    figure
    
    if class(spiketimes) == "cell"
        hold on
        for i = 1:length(spiketimes)
            indiv_trial = spiketimes{i};
            sz = size(indiv_trial);
            plot(spiketimes{i}, i*ones(sz), 'k.', 'MarkerSize', 2)
        end
        hold off
        
    elseif class(spiketimes) == "double" %i.e. if single trial
        sz = size(spiketimes);
        plot(spiketimes, ones(sz), 'k.', 'MarkerSize', 2)
    end   
    
    xlabel('time /ms')
    ylabel('trial')
    title('Raster plot of spike times')
end

