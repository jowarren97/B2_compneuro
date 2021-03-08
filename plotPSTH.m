function [] = plotPSTH(spiketimes, dt)
%UNTITLED2 Summary of this function goes here
%   t_end, dt in milliseconds

    if class(spiketimes) == "cell"
        allspiketimes = [spiketimes{:}];
        t_end = ceil(max(allspiketimes)/100)*100; %round t ax up to nearest 100
        n_trials = length(spiketimes);
    elseif class(spiketimes) == "double"
        t_end = ceil(max(spiketimes)/100)*100;
        n_trials = 1;
        allspiketimes = spiketimes;
    end
    
    figure('Position', [100, 100, 600, 400])
    edges = 0:dt:t_end;
    response = histc(allspiketimes, edges);
    response = response*1000/dt; %normalise for time
    response = response/n_trials; %normalise for number of trials
    bar(edges, response, 'histc');
    xlabel('time /ms')
    ylabel('spike rate /s^{-1}')
    xlim([0, t_end])
    title('PSTH')
    set(gca, 'Fontsize', 14)
end

