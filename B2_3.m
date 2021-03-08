%%
clear all;
load data/spiketimes-revcorr.1.mat

%%
last_sample_time = 6.25*(size(stimulus_ft, 2)-1);
t = 0:6.25:last_sample_time;
f = logspace(log10(500), log10(22400), 34);
%% Q4
imagesc(t,f, stimulus_ft)
colormap(jet)
axis xy
xlabel('time /ms')
ylabel('frequency /Hz')
cb = colorbar;
title('Stimulus spectrogram')
cb.Title.String = 'Amplitude /dB';
set(gca, 'Fontsize', 14)
saveas(gcf, 'figs/png/B2_3_q4.png')
saveas(gcf, 'figs/mat/B2_3_q4.fig')
%% Q6
single_trial = spiketimes{1};
plotraster(single_trial)
saveas(gcf, 'figs/png/B2_3_q6.png')
saveas(gcf, 'figs/mat/B2_3_q6.fig')

%% Q8
plotraster(spiketimes)
saveas(gcf, 'figs/png/B2_3_q8.png')
saveas(gcf, 'figs/mat/B2_3_q8.fig')

%% Q9
plotPSTH(spiketimes, 6.25)
saveas(gcf, 'figs/png/B2_3_q9.png')
saveas(gcf, 'figs/mat/B2_3_q9.fig')

%% Q11
computeSTA(spiketimes, stimulus_ft, 15)
saveas(gcf, 'figs/png/B2_3_q11.png')
saveas(gcf, 'figs/mat/B2_3_q11.fig')

%% Q12
allspiketimes = [spiketimes{:}]
kernel= separablekernel(stimulus_ft, histc(allspiketimes, t), 15);
figure
imagesc(kernel)
saveas(gcf, 'figs/png/B2_3_q12.png')
saveas(gcf, 'figs/mat/B2_3_q12.fig')

