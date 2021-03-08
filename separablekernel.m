function kernel = separablekernel(X_ft, y_t, n_h)
% function kernel = separablekernel(X_ft, y_t, n_h)
% 
% Calculate separable kernel.
%
% Inputs:
%  X_ft -- stimulus matrix, frequency x time
%  y_t -- neuronal response vector
%  n_h -- number of history steps required
% 
% Outputs:
%  kernel -- the separable kernel

y_t = y_t(:);

n_f = size(X_ft, 1);
n_t = size(X_ft, 2);

% pad with zeros
X_ft = [zeros(n_f, n_h) X_ft];

n_t_pad = size(X_ft, 2);

% preallocate
X_fht = zeros(n_f, n_h, n_t_pad);

for ii = 1:n_h
  X_fht(:,ii,:) = reshape(circshift(X_ft, [0 n_h-ii]), [n_f, 1, n_t_pad]);
end

X_fht = X_fht(:, :, n_h+1:end);

% constant terms
X_fht(end+1, end+1, :) = 1;

[n_f, n_h, n_t] = size(X_fht);

k_f = ones(n_f, 1);
k_h = ones(n_h, 1);

for ii = 1:15
 yh = X_fht.*repmat(k_f, [1 n_h n_t]);
 yh = squeeze(sum(yh, 1));
 k_h = yh'\y_t;

 yh = X_fht.*repmat(k_h', [n_f 1 n_t]);
 yh = squeeze(sum(yh, 2));
 k_f = yh'\y_t;
end

% separate out constant terms
kernel.c_f = k_f(end);
kernel.k_f = k_f(1:end-1);
kernel.c_h = k_h(end);
kernel.k_h = k_h(1:end-1);

kernel = kernel.k_f * kernel.k_h';
