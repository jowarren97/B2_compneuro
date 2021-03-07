function [] = plotPSTH(spiketimes, t_end, dt)
%UNTITLED2 Summary of this function goes here
%   t_end, dt in milliseconds
    if class(spiketimes) == "cell"
        n_trials = length(spiketimes);
        allspiketimes = [spiketimes{:}];
    elseif class(spiketimes) == "double"
        n_trials = 1;
        allspiketimes = spiketimes;
    end
        
    edges = 0:dt:t_end;
    response = histc(allspiketimes, edges);
    response = response*1000/dt; %normalise for time
    response = response/n_trials; %normalise for number of trials
    bar(edges, response, 'histc');
    xlabel('time/ ms')
    ylabel('spike rate/ 1/s')
    xlim([0, t_end])
    title('PSTH from')
end

