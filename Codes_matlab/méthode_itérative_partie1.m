clear; close all; clc;

%% ============================================================
%  Script unique : methode 2 points iterative (AM Alg1)
%  Figures approx du papier :
%   - RMSE vs delta  (proche Fig. 5)
%   - RMSE vs SNR    (proche Fig. 7)
%% ============================================================

rng(1);

%% ============================================================
%  FIGURE A : RMSE vs delta  (approx Fig. 5)
%  Papier : N = 512, SNR = 0 dB
%% ============================================================
N      = 512;
fs     = 512e3;
A      = 1;
SNR_dB = 0;
Nruns  = 1000;

delta_vals = -0.49:0.02:0.49;

rmse_am_iter1 = zeros(size(delta_vals));
rmse_am_iter2 = zeros(size(delta_vals));
crlb_rmse     = zeros(size(delta_vals));

n = 0:N-1;

for id = 1:length(delta_vals)
    delta_true = delta_vals(id);
    f0 = (N/4 + delta_true) * fs / N;   % comme indique dans le papier

    est1 = zeros(1, Nruns);
    est2 = zeros(1, Nruns);

    for r = 1:Nruns
        theta0 = 2*pi*rand;

        % signal complexe, comme dans le modele du papier
        s = A * exp(1j*(2*pi*f0*n/fs + theta0));

        % bruit complexe blanc, SNR = A^2 / sigma^2
        sigma2 = A^2 / (10^(SNR_dB/10));
        z = sqrt(sigma2/2) * (randn(1,N) + 1j*randn(1,N));

        x = s + z;

        % 1 iteration
        delta_hat_1 = am_alg1_estimator(x, 1);

        % 2 iterations
        delta_hat_2 = am_alg1_estimator(x, 2);

        est1(r) = delta_hat_1;
        est2(r) = delta_hat_2;
    end

    rmse_am_iter1(id) = sqrt(mean((est1 - delta_true).^2));
    rmse_am_iter2(id) = sqrt(mean((est2 - delta_true).^2));

    % CRLB sur delta, d'apres Eq. (23) du PDF
    % Var(delta) ~= 6/( (2*pi)^2 ) * sigma^2/A^2 * 1/N
    sigma2 = A^2 / (10^(SNR_dB/10));
    var_crlb_delta = 6 * sigma2 / (A^2 * (2*pi)^2 * N);
    crlb_rmse(id)  = sqrt(var_crlb_delta);
end

figure;
plot(delta_vals, rmse_am_iter1, 'o', 'MarkerSize', 4); hold on;
plot(delta_vals, rmse_am_iter2, '-', 'LineWidth', 1.6);
plot(delta_vals, crlb_rmse, '--', 'LineWidth', 1.6);
grid on;
xlabel('\delta');
ylabel('RMSE');
title('AM Alg1 (2 points itératif) : RMSE vs \delta');
legend('AM itération 1', 'AM itération 2', 'CRLB approx', 'Location', 'northwest');

%% ============================================================
%  FIGURE B : RMSE vs SNR  (approx Fig. 7)
%  Papier : N = 16, delta = 0.2
%% ============================================================
N      = 16;
fs     = 16e3;
A      = 1;
delta_true = 0.2;
Nruns  = 5000;
SNR_vals = 0:5:90;

rmse_iter2 = zeros(size(SNR_vals));
crlb_snr   = zeros(size(SNR_vals));

n = 0:N-1;
f0 = (N/4 + delta_true) * fs / N;

for is = 1:length(SNR_vals)
    SNR_dB = SNR_vals(is);
    est2 = zeros(1, Nruns);

    for r = 1:Nruns
        theta0 = 2*pi*rand;

        s = A * exp(1j*(2*pi*f0*n/fs + theta0));

        sigma2 = A^2 / (10^(SNR_dB/10));
        z = sqrt(sigma2/2) * (randn(1,N) + 1j*randn(1,N));

        x = s + z;

        est2(r) = am_alg1_estimator(x, 2);
    end

    rmse_iter2(is) = sqrt(mean((est2 - delta_true).^2));

    var_crlb_delta = 6 * sigma2 / (A^2 * (2*pi)^2 * N);
    crlb_snr(is)   = sqrt(var_crlb_delta);
end

figure;
semilogy(SNR_vals, rmse_iter2, 'o-', 'LineWidth', 1.6, 'MarkerSize', 4); hold on;
semilogy(SNR_vals, crlb_snr, '-', 'LineWidth', 1.6);
grid on;
xlabel('SNR (dB)');
ylabel('RMSE');
title('AM Alg1 (2 points itératif, 2 itérations) : RMSE vs SNR');
legend('AM itération 2', 'CRLB approx', 'Location', 'southwest');

%% ============================================================
%  FONCTIONS LOCALES
%% ============================================================

function delta_hat = am_alg1_estimator(x, Q)
    % AM Alg1 du tableau 3 du papier
    % Initialisation : delta_hat = 0
    % Puis Q iterations
    %
    % x : signal complexe
    % Q : nombre d'iterations

    N = length(x);
    n = 0:N-1;

    % FFT pour trouver m (coarse estimation)
    Xfft = fft(x);
    Xh = abs(Xfft(1:floor(N/2)));
    [~, kmax] = max(Xh);

    % en indice MATLAB, frequence coarse = kmax-1 bins
    m_hat = kmax - 1;

    delta_hat = 0;

    for q = 1:Q
        p1 = delta_hat + 0.5;
        p2 = delta_hat - 0.5;

        Xp1 = sum(x .* exp(-1j*2*pi*n*(m_hat + p1)/N));
        Xp2 = sum(x .* exp(-1j*2*pi*n*(m_hat + p2)/N));

        denom = (Xp1 - Xp2);
        if abs(denom) < 1e-14
            break;
        end

        delta_corr = 0.5 * real((Xp1 + Xp2) / denom);
        delta_hat = delta_hat + delta_corr;

        % sécurité
        if abs(delta_corr) < 1e-10
            break;
        end
    end
end