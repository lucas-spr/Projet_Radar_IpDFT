clear; close all; clc;

%% ============================================================
%  RMSE/CRLB vs delta pour la methode iterative 2 points (AM Alg1)
%  Inspiration directe de la Fig. 5 du papier
%  Parametres du papier : N = 512, SNR = 0 dB
%% ============================================================

rng(1);

% Parametres
N      = 512;
fs     = 512e3;
A      = 1;
SNR_dB = 0;
Nruns  = 1000;

delta_vals = -0.49:0.02:0.49;
n = 0:N-1;

rmse_iter1 = zeros(size(delta_vals));
rmse_iter2 = zeros(size(delta_vals));
crlb_rmse  = zeros(size(delta_vals));

% bruit complexe du modele du papier : SNR = A^2 / sigma^2
sigma2 = A^2 / (10^(SNR_dB/10));

% CRLB sur delta (approximation utilisee dans le papier)
% Var(delta_hat)_CR ~= 6 / ((2*pi)^2 * N * SNR)
var_crlb_delta = 6 * sigma2 / (A^2 * (2*pi)^2 * N);
rmse_crlb_abs  = sqrt(var_crlb_delta);

for id = 1:length(delta_vals)

    delta_true = delta_vals(id);

    % comme indique dans le papier :
    % f0 = (N/4 + delta)*fs/N
    f0 = (N/4 + delta_true) * fs / N;

    est1 = zeros(1, Nruns);
    est2 = zeros(1, Nruns);

    for r = 1:Nruns

        theta0 = 2*pi*rand;

        % signal complexe du papier
        s = A * exp(1j*(2*pi*f0*n/fs + theta0));

        % bruit complexe blanc gaussien
        z = sqrt(sigma2/2) * (randn(1,N) + 1j*randn(1,N));

        x = s + z;

        % estimateur AM itération 1
        est1(r) = am_alg1_estimator(x, 1);

        % estimateur AM itération 2
        est2(r) = am_alg1_estimator(x, 2);
    end

    rmse_iter1(id) = sqrt(mean((est1 - delta_true).^2));
    rmse_iter2(id) = sqrt(mean((est2 - delta_true).^2));
    crlb_rmse(id)  = rmse_crlb_abs;
end

% Normalisation RMSE / CRLB
ratio_iter1 = rmse_iter1 ./ crlb_rmse;
ratio_iter2 = rmse_iter2 ./ crlb_rmse;
ratio_crlb  = crlb_rmse  ./ crlb_rmse;

%% ============================================================
%  TRACE
%% ============================================================
figure;
plot(delta_vals, ratio_iter1, 'ks', 'MarkerSize', 5, 'LineWidth', 1.0); hold on;
plot(delta_vals, ratio_iter2, 'ko', 'MarkerSize', 5, 'LineWidth', 1.0);
plot(delta_vals, ratio_crlb,  'k-', 'LineWidth', 1.5);

grid on;
xlabel('\delta');
ylabel('RMSE / CRLB');
title('RMSE/CRLB en fonction de \delta - AM estimator (2 points itératif)');
legend('AM estimator (iteration 1, simulated)', ...
       'AM estimator (iteration 2, simulated)', ...
       'CRLB', ...
       'Location', 'north');

xlim([-0.5 0.5]);
ylim([0.95 1.9]);

%% ============================================================
%  FONCTION LOCALE : AM Alg1
%  Tableau 3 du papier
%% ============================================================
function delta_hat = am_alg1_estimator(x, Q)

    N = length(x);
    n = 0:N-1;

    % recherche du pic principal
    Xfft = fft(x);
    Xh = abs(Xfft(1:floor(N/2)));
    [~, kmax] = max(Xh);

    % index en bins
    m_hat = kmax - 1;

    % initialisation
    delta_hat = 0;

    for q = 1:Q

        % evaluation de la DFT a m + delta_hat +/- 0.5
        p_plus  = delta_hat + 0.5;
        p_minus = delta_hat - 0.5;

        Xp = sum(x .* exp(-1j*2*pi*n*(m_hat + p_plus )/N));
        Xm = sum(x .* exp(-1j*2*pi*n*(m_hat + p_minus)/N));

        denom = (Xp - Xm);
        if abs(denom) < 1e-14
            break;
        end

        % formule AM Alg1
        delta_corr = 0.5 * real((Xp + Xm) / denom);

        delta_hat = delta_hat + delta_corr;

        % critere d'arret
        if abs(delta_corr) < 1e-12
            break;
        end
    end
end