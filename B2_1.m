%% PARAMS

clear all
tau_m = 10; %ms
V_rest = -70; %mV
V_reset = -70
V_thresh = -40; %mV
R_m = 10; %MOhm
dt = 1; %ms
T = 1000; %ms
G = 1/R_m;

t = 0:dt:T;
steps = T/dt;

%% QU 1
I = 3.1; %nA
V = V_rest*ones(1,steps+1);

for j = 1:steps
    V(j+1) = V(j) + (1/tau_m)*dt*(R_m*I - (V(j)-V_rest));
    if V(j+1) >= V_thresh
        V(j+1) = V_reset;
    end        
end

figure

hold on
plot(t, V, 'lineWidth', 1.5)
plot(t, V_rest*ones(size(t)), 'k--', 'lineWidth', 1.5)
plot(t, V_thresh*ones(size(t)), 'r--', 'lineWidth', 1.5)
hold off
legend('Membrane potential', 'Resting potential', 'Threshold potential')

xlim([0,1000])
ylim([V_rest-10, V_thresh+10])
set(gca, 'Fontsize', 14)
xlabel('time /ms')
ylabel('Membrane Voltage /mV')

title(sprintf('LIF simulation, input current = %.1fnA', I))

saveas(gcf, 'figs/png/B2_1_q1.png')
saveas(gcf, 'figs/mat/B2_1_q1.fig')

%% QU 3
I = 2.9; %nA
V = V_rest*ones(1,steps+1);

for j = 1:steps
    V(j+1) = V(j) + (1/tau_m)*dt*(R_m*I - (V(j)-V_rest));
    if V(j+1) >= V_thresh
        V(j+1) = V_reset;
    end        
end

figure

hold on
plot(t, V, 'lineWidth', 1.5)
plot(t, V_rest*ones(size(t)), 'k--', 'lineWidth', 1.5)
plot(t, V_thresh*ones(size(t)), 'r--', 'lineWidth', 1.5)
hold off
legend('Membrane potential', 'Resting potential', 'Threshold potential')

xlim([0,1000])
ylim([V_rest-10, V_thresh+10])
set(gca, 'Fontsize', 14)
xlabel('time /ms')
ylabel('Membrane Voltage /mV')

title(sprintf('LIF simulation, input current = %.1fnA', I))

%legend(sprintf('I = %.1fnA', I))
saveas(gcf, 'figs/png/B2_1_q3.png')
saveas(gcf, 'figs/mat/B2_1_q3.fig')

%%
I_trials = 2:0.1:5;
n_trials = length(I_trials);
V_trials = V_rest*ones(n_trials, steps+1); % nb trials x duration
spike_counts = zeros(1, n_trials);

for k = 1:n_trials
    I = I_trials(k);
    for j = 1:steps
        V_trials(k,j+1) = V_trials(k,j) + (1/tau_m)*dt*(R_m*I - (V_trials(k,j)-V_rest));
        if V_trials(k,j+1) >= V_thresh
            V_trials(k,j+1) = V_reset;
            spike_counts(k) = spike_counts(k) + 1;
        end            
    end
end

rates = spike_counts*1000/T;

figure
bar(I_trials, rates)
set(gca, 'Fontsize', 14)
title('Firing rate of LIF model as function of input current')
xlabel('Input current /nA')
ylabel('Firing rate /s^{-1}')
%legend(sprintf('I = %.1fnA', I))
saveas(gcf, 'figs/png/B2_1_q4.png')
saveas(gcf, 'figs/mat/B2_1_q4.fig')

%%
clear all
%neuron params
tau_m = 20; %ms
V_rest = -70; %mV
V_reset = -80;
V_thresh = -54; %mV
R_m_I_e = 18; %mV
dt = 1; %ms
T = 1000; %ms
%syn params
R_m_G_s = 0.15;
P_max = 0.5;
tau_s = 10; %ms

t = 0:dt:T;
steps = T/dt;

n_neurons = 2;
%%
E_s = 0;

V_init = V_reset+(V_thresh-V_reset)*rand([n_neurons,1]);
V = repmat(V_init, 1,steps+1);

%connectivity = [0, 1; 1, 0]; % [1->1 1->2; 2->1 2->2]
connectivity = ones(n_neurons, n_neurons) .* ~eye(n_neurons, n_neurons) %remore self connections
%P = zeros([2,steps+1]); %P = [P_1->2; P_2->1]
P = zeros([n_neurons, steps+1]);
syn_inp = zeros([n_neurons, steps+1]);

for j = 1:steps
    syn_mat = connectivity.*P(:,j);
    syn_inp(:,j) = - R_m_G_s*syn_mat'*(V(:,j)-E_s);
    V(:,j+1) = V(:,j) + (1/tau_m)*dt*(R_m_I_e - (V(:,j)-V_rest) + syn_inp(:,j));
    P(:,j+1) = P(:,j)*(1 - (1/tau_s)*dt);
    spiked_bool = V(:,j+1) >= V_thresh
    V(spiked_bool, j+1) = V_reset;
    P(spiked_bool, j+1) = P_max;
    spiketrains(:,j+1) = spiked_bool;
end

fig=figure('Position', [600,600,1000,600])

subplot(2,1,1)
hold on
plot(t, V, 'lineWidth', 1.5)
plot(t, V_rest*ones(size(t)), 'k--', 'lineWidth', 1.5)
plot(t, V_thresh*ones(size(t)), 'm--', 'lineWidth', 1.5)
plot(t, V_reset*ones(size(t)), 'b--', 'lineWidth', 1.5)
hold off
legend('Neuron 1', 'Neuron 2', 'Resting potential', 'Threshold potential', 'Reset potential')

xlim([0,1000])
ylim([-85, -50])
set(gca, 'Fontsize', 14)
xlabel('time /ms')
ylabel('Membrane Voltage /mV')

subplot(2,1,2)
plot(t, syn_inp, 'lineWidth', 1.5)
legend('Neuron 1', 'Neuron 2')
xlim([0,1000])
set(gca, 'Fontsize', 14)
xlabel('time /ms')
ylabel('Synaptic current /nA')

ax=axes(fig,'visible','off'); 
ax.Title.Visible='on';
title(ax, sprintf("Simulation of excitatory LIF population of %d neurons", n_neurons))
set(gca, 'Fontsize', 14)

%legend(sprintf('I = %.1fnA', I))
saveas(gcf, 'figs/png/B2_1_q5a.png')
saveas(gcf, 'figs/mat/B2_1_q5a.fig')

%%
E_s = -80;

V_init = V_reset+(V_thresh-V_reset)*rand([n_neurons,1]);
V = repmat(V_init, 1,steps+1);

%connectivity = [0, 1; 1, 0]; % [1->1 1->2; 2->1 2->2]
connectivity = ones(n_neurons, n_neurons) .* ~eye(n_neurons, n_neurons) %remore self connections
%P = zeros([2,steps+1]); %P = [P_1->2; P_2->1]
P = zeros([n_neurons, steps+1]);
syn_inp = zeros([n_neurons, steps+1]);

for j = 1:steps
    syn_mat = connectivity.*P(:,j)
    syn_inp(:,j) = - R_m_G_s*syn_mat'*(V(:,j)-E_s);
    V(:,j+1) = V(:,j) + (1/tau_m)*dt*(R_m_I_e - (V(:,j)-V_rest) + syn_inp(:,j));
    P(:,j+1) = P(:,j)*(1 - (1/tau_s)*dt);
    spiked_bool = V(:,j+1) >= V_thresh
    V(spiked_bool, j+1) = V_reset;
    P(spiked_bool, j+1) = P_max;
    spiketrains(:,j+1) = spiked_bool;
end

fig=figure('Position', [600,600,1000,600])

subplot(2,1,1)
hold on
plot(t, V, 'lineWidth', 1.5)
plot(t, V_rest*ones(size(t)), 'k--', 'lineWidth', 1.5)
plot(t, V_thresh*ones(size(t)), 'm--', 'lineWidth', 1.5)
plot(t, V_reset*ones(size(t)), 'b--', 'lineWidth', 1.5)
hold off
legend('Neuron 1', 'Neuron 2', 'Resting potential', 'Threshold potential', 'Reset potential')

xlim([0,1000])
ylim([-85, -50])
set(gca, 'Fontsize', 14)
xlabel('time /ms')
ylabel('Membrane Voltage /mV')

subplot(2,1,2)
plot(t, syn_inp, 'lineWidth', 1.5)
legend('Neuron 1', 'Neuron 2')
xlim([0,1000])
set(gca, 'Fontsize', 14)
xlabel('time /ms')
ylabel('Synaptic current /nA')

ax=axes(fig,'visible','off'); 
ax.Title.Visible='on';
title(ax, sprintf("Simulation of inhibitory LIF population of %d neurons", n_neurons))
set(gca, 'Fontsize', 14)

%legend(sprintf('I = %.1fnA', I))
saveas(gcf, 'figs/png/B2_1_q5b.png')
saveas(gcf, 'figs/mat/B2_1_q5b.fig')


%%
clear all
close all

%net params
learning = false;
n_neurons = 100; %nb recurr conn neurons
n_inp = 1; %nb input channels
conn_rec = ones([n_neurons, n_neurons]);
conn_fwd = ones([n_inp, n_neurons]);
regL2 = 0;
regL1 = 0;

%neuron params
tau_m = 20; %ms
V_reset = -0.5;
V_thresh = 0.5;
sigma_V = 10^-2;

%sim params
dt = 1; %ms
T = 4000; %ms
t = 0:dt:T;
steps = T/dt;

%initialisation
V_init = V_reset+(V_thresh-V_reset)*rand([n_neurons, 1]);
V = repmat(V_init, 1, steps+1);

spiketrains = zeros([n_neurons, steps+1]);
rates = zeros([n_neurons, steps+1]); %filtered/low-pass of spiketrains
output = zeros([n_inp, steps+1]);

W_fwd = 0.5*random('normal', 0, 1, [n_inp, n_neurons]);
W_rec = - W_fwd' * W_fwd - regL2*(1/tau_m)*eye([n_neurons, n_neurons]);

V_thresh = 0.5 * diag(W_fwd' * W_fwd) + regL1*(1/tau_m) + regL2*(1/tau_m)^2;

figure
imagesc(W_fwd)
title('Initial feedforward weight matrix')

figure
imagesc(W_rec)
title('Initial recurrent weight matrix')

omega = 1/4000;
amplitude = 2;
x = amplitude * sin(2*pi*omega*t);
x_dot = amplitude * 2*pi*omega*cos(2*pi*omega*t);
c = x + (1/tau_m)*x_dot;

for j = 1:steps
    %c(:,j) = x(:,j+1) - x(:,j) + (1/tau_m)*dt*x(:,j); 
    
    syn_mat_fwd = conn_fwd.*W_fwd;
    syn_mat_rec = conn_rec.*W_rec;
    
    syn_inp(:,j) = dt*(1/tau_m) * syn_mat_fwd' * c(:,j) + syn_mat_rec * spiketrains(:,j);
    
    V_noise = sigma_V * random('normal', 0, 1, [n_neurons,1]); 
    
    V(:,j+1) = V(:,j) - (1/tau_m)*dt*V(:,j) + syn_inp(:,j) + V_noise;

    [V_above_thresh, idx_spiked] = max(V(:,j+1) - V_thresh); %only allow 1 neuron to spike per timestep (compensate for discrete simulation)
                                        %neuron firing should inhibit neighbours
                                        %want to take neuron whose V is most above threshold and assume spiked first
    %spiked_bool = V(:,j+1) >= V_thresh;
    if V_above_thresh >= 0
        spiketrains(idx_spiked,j+1) = 1;
        
        if learning == true %not sure if works
            dw = -2*(V(:,j) + regL2*rate(:,j)) - W_rec(:,idx_spiked);
            dw(idx_spiked) = dw(idx_spiked) - regL2;
            W_rec(:,idx_spiked) = W_rec(:,idx_spiked) + lr_rec*dw;
        end

    end
    rates(:,j+1) = (1 - dt/tau_m) * rates(:,j) + spiketrains(:,j);
    output(:,j+1) = (1-dt/tau_m) * output(:,j) + W_fwd * spiketrains(:,j+1);

end


figure
plot(t, V);
title('Membrane voltage of neurons')
% figure
% imagesc(W_fwd(:,:,steps+1));
% figure
% imagesc(W_rec(:,:,steps+1));
figure
hold on
plot(t, output)
plot(t, x, '--', 'linewidth', 3)
title(sprintf('Output of network of %d neurons', n_neurons))
legend('Network output', 'Input signal')
xlabel('time /ms')
set(gca, 'Fontsize', 14)

saveas(gcf, 'figs/png/B2_1_q6.png')

