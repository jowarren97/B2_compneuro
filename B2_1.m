%% QU 1

tau_m = 10; %ms
V_rest = -70; %mV
V_thresh = -40; %mV
R_m = 10; %MOhm
dt = 1; %ms
T = 1000; %ms
G = 1/R_m;

t = 0:dt:T;
steps = T/dt;

%%
I = 3.1; %nA
V = V_rest*ones(1,steps+1);

for j = 1:steps
    V(j+1) = V(j) + (1/tau_m)*dt*(I - G*(V(j)-V_rest));
    if V(j+1) >= V_thresh
        V(j+1) = V_rest;
    end        
end

figure()
plot(V)

%%
I = 2.9; %nA
V = V_rest*ones(1,steps+1);

for j = 1:steps
    V(j+1) = V(j) + (1/tau_m)*dt*(I - G*(V(j)-V_rest));
    if V(j+1) >= V_thresh
        V(j+1) = V_rest;
    end        
end

figure()
plot(V)

%%
I_trials = 2:0.1:5;
n_trials = length(I_trials);
V_trials = V_rest*ones(n_trials, steps+1); % nb trials x duration
spike_counts = zeros(1, n_trials);

for k = 1:n_trials
    I = I_trials(k);
    for j = 1:steps
        V_trials(k,j+1) = V_trials(k,j) + (1/tau_m)*dt*(I - G*(V_trials(k,j)-V_rest));
        if V_trials(k,j+1) >= V_thresh
            V_trials(k,j+1) = V_rest;
            spike_counts(k) = spike_counts(k) + 1;
        end            
    end
end

rates = spike_counts/T;
%%
figure()
plot(I_trials, rates)