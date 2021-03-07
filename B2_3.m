%%
load spiketimes-revcorr.1.mat

%%
last_sample_time = 6.25*(size(stimulus_ft, 2)-1);
t = 0:6.25:last_sample_time;
f = logspace(log10(500), log10(22400), 34);
%%
imagesc(t,f, stimulus_ft)
colormap(jet)
axis xy
xlabel('time /ms')
ylabel('frequency /Hz')
cb = colorbar;
title('Stimulus spectrogram')
cb.Title.String = 'Amplitude /dB'

%%
single_trial = spiketimes{1};
plotraster(single_trial)

%%
plotraster(spiketimes)

%%
plotPSTH(spiketimes, last_sample_time, 6.25)

%%
computeSTA(spiketimes, stimulus_ft, 15)
