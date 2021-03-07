load spiketimes.mat

%% 
sz = size(single_trial);

first_spike = single_trial(1);
first_ten_spikes = single_trial(1:10);
plot(single_trial, ones(sz), 'o')
xlabel('time /ms')
title('Raster plot of single trial')

%% 
edges = 0:10:1000;
response = histc(single_trial, edges);
response = response/0.01;
bar(edges, response, 'histc');
xlabel('time/ ms')
ylabel('spike count')
title('PSTH from a single trial')

%% 
first_hist = condition1(1,:);
figure()
bar(edges, first_hist);

figure()
for ii = 1:20
 subplot(5,4,ii);
 bar(edges, condition1(ii,:)/0.01,'histc');
end

%% 
figure()
for ii = 1:length(condition1)
 subplot(5,4,ii);
 bar(edges, condition1(:,ii)/0.01,'histc');
end
%??????

%% 
mn = zeros(1,101);
for ii = 1:20
    mn = mn + condition1(ii,:)/0.01;
end
mn = mn/20;
bar(edges, mn, 'histc')

%% 
figure()
mu1 = mean(condition1)/0.01;
sigma1 = std(condition1)/0.01;
standard_error1 = sigma/sqrt(20);
bar(edges, mu1, 'histc')
hold all
errorbar(edges+5, mu1, standard_error1, 'LineStyle', 'none')
hold off
ylim([0,800])

figure()
mu2 = mean(condition2)/0.01;
sigma2 = std(condition2)/0.01;
standard_error2 = sigma/sqrt(20);
bar(edges, mu2, 'histc')
hold all
errorbar(edges+5, mu2, standard_error2, 'LineStyle', 'none')
hold off
ylim([0,800])

figure()
hold all
errorbar(edges+5, mu1, standard_error1)
errorbar(edges+5, mu2, standard_error2)
hold off






