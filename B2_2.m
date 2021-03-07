load data/spiketimes.mat

%% Q3
sz = size(single_trial);

first_spike = single_trial(1);
first_ten_spikes = single_trial(1:10);
figure('Position', [100, 100, 600, 200])
plot(single_trial, ones(sz), 'o')
xlabel('time /ms')
title('Raster plot of single trial')
set(gca, 'Fontsize', 14, 'ytick',[])
saveas(gcf, 'figs/png/B2_2_q3.png')
saveas(gcf, 'figs/mat/B2_2_q3.fig')

%% Q4
edges = 0:10:1000;
response = histc(single_trial, edges);
response = response/0.01;
figure;
bar(edges, response, 'histc');
xlabel('time /ms')
ylabel('Spike rate /s^{-1}')
title('PSTH from a single trial')
set(gca, 'Fontsize', 14)
saveas(gcf, 'figs/png/B2_2_q4.png')
saveas(gcf, 'figs/mat/B2_2_q4.fig')

%% Q7
first_hist = condition1(1,:);
figure;
bar(edges, first_hist);

fig = figure;
for ii = 1:20
 subplot(5,4,ii);
 bar(edges, condition1(ii,:)/0.01,'histc');
 ylim([0, 1000]);
 if ii <= 16
     set(gca, 'xtick', [])
 end
 if mod(ii,4) ~= 1
     set(gca, 'ytick', [])
 end
end
ax=axes(fig,'visible','off'); 
ax.Title.Visible='on';
ax.XLabel.Visible='on';
ax.YLabel.Visible='on';
ylabel(ax,'Spike rate /s^{-1}');
xlabel(ax,'time /ms');
title(ax,'PSTH plots for 20 trials');
set(gca, 'Fontsize', 14)
saveas(gcf, 'figs/png/B2_2_q7.png')
saveas(gcf, 'figs/mat/B2_2_q7.fig')

%% Q8
fig=figure('Position', [10 10 800 600]);
for ii = 1:length(condition1)-1
 subplot(10,10,ii);
 bar(1:20, condition1(:,ii)/0.01,'histc');
 ylim([0, 1000]);
 title(sprintf('%d-%dms', (ii-1)/0.1, ii/0.1))
 if ii <= 90
     set(gca, 'xtick', [])
 end
 if mod(ii,10) ~= 1
     set(gca, 'ytick', [])
 end
end
ax=axes(fig,'visible','off'); 
ax.Title.Visible='on';
ax.XLabel.Visible='on';
ax.YLabel.Visible='on';
ylabel(ax,'Spike rate /s^{-1}');
xlabel(ax,'Trial number');
title(ax,'Histograms of activity across trials');
set(gca, 'Fontsize', 14)
saveas(gcf, 'figs/png/B2_2_q8.png')
saveas(gcf, 'figs/mat/B2_2_q8.fig')

%% 
mn = zeros(1,101);
for ii = 1:20
    mn = mn + condition1(ii,:)/0.01;
end
mn = mn/20;
bar(edges, mn, 'histc')

xlabel('time /ms')
ylabel('Spike rate /s^{-1}')
title('PSTH averaged across 20 trials')
set(gca, 'Fontsize', 14)
saveas(gcf, 'figs/png/B2_2_q9.png')
saveas(gcf, 'figs/mat/B2_2_q9.fig')

%% 
% figure()
% mu1 = mean(condition1)/0.01;
% sigma1 = std(condition1)/0.01;
% standard_error1 = sigma1/sqrt(20);
% bar(edges, mu1, 'histc')
% hold all
% errorbar(edges+5, mu1, standard_error1, 'LineStyle', 'none')
% hold off
% ylim([0,800])
% 
% figure()
% mu2 = mean(condition2)/0.01;
% sigma2 = std(condition2)/0.01;
% standard_error2 = sigma2/sqrt(20);
% bar(edges, mu2, 'histc')
% hold all
% errorbar(edges+5, mu2, standard_error2, 'LineStyle', 'none')
% hold off
% ylim([0,800])

figure('Position', [10 10 800 600])
hold all
errorbar(edges+5, mu1, standard_error1, 'lineWidth', 1)
errorbar(edges+5, mu2, standard_error2, 'lineWidth', 1)
hold off

xlim([0 1000])
xlabel('time /ms')
ylabel('Spike rate /s^{-1}')
set(gca, 'Fontsize', 14)
legend('Condition 1', 'Condition 2')
set(gca, 'Fontsize', 14)
title('PSTH plots for conditions 1 and 2')
saveas(gcf, 'figs/png/B2_2_q10.png')
saveas(gcf, 'figs/mat/B2_2_q10.fig')

%%
col1 = ismember(sigma1, max(sigma1(:)));
col1 = condition1(:,col1);
col2 = ismember(sigma2, max(sigma2(:)));
col2 = condition2(:,col2);

%????

%%
spiketimes1{1}
plotraster(spiketimes1)
set(gca, 'Fontsize', 14)
saveas(gcf, 'figs/png/B2_2_q13.png')
saveas(gcf, 'figs/mat/B2_2_q13.fig')

%%
plotPSTH(spiketimes1, 10)
set(gca, 'Fontsize', 14)
title("PSTH of 'spiketimes1'")
saveas(gcf, 'figs/png/B2_2_q14.png')
saveas(gcf, 'figs/mat/B2_2_q14.fig')

