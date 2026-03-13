clear; clc; close all;
fprintf('=== SIMULATION SUR MULTIPLES DELTAS ===\n');

% 1. Paramètres globaux
M = 1024;               
fs = 1000;              
t = 0:1/fs:(M-1)/fs;    
Ak = 2;                 
bin_k = 13;             
H = 2;                  % Fenêtre de Hann (H=2)

SNR_dB_range = 40:10:100; 
nb_runs = 1000;         % Itérations Monte-Carlo

deltas_test = [0.1, -0.4, 0.3, -0.25];

% Création de la figure (Format 2x2)
figure('Name', 'Variances pour différents \delta_k', 'Color', 'white', 'Position', [100, 100, 1000, 800]);

%% 2. Boucle sur les différents Deltas
for d_idx = 1:length(deltas_test)
    delta_k = deltas_test(d_idx);
    fk = (bin_k + delta_k) * (fs / M);
    
    fprintf('Traitement en cours pour delta_k = %.2f...\n', delta_k);
    
    var_crlb_th = zeros(1, length(SNR_dB_range));
    var_ipdft_th = zeros(1, length(SNR_dB_range));
    var_ipdft_sim = zeros(1, length(SNR_dB_range));
    
    % Boucle sur les SNR
    for idx_snr = 1:length(SNR_dB_range)
        SNR_dB = SNR_dB_range(idx_snr);
        
        % Puissance et Bruit
        Puissance_Signal = (Ak^2) / 2;
        SNR_lineaire = 10^(SNR_dB / 10);
        sigma_n_carre = Puissance_Signal / SNR_lineaire;
        sigma_n = sqrt(sigma_n_carre);
        
        % THÉORIE
        var_crlb_th(idx_snr) = (6 / pi^2) * (sigma_n_carre / (M * Ak^2));
        var_ipdft_th(idx_snr) = calculer_variance_eq13(delta_k, Ak, M, sigma_n, H);
        
        % SIMULATION PRATIQUE
        deltas_estimes = zeros(1, nb_runs);
        for run = 1:nb_runs
            phi = rand * 2 * pi;
            sig_pur = Ak * cos(2*pi*fk*t + phi);
            sig_bruite = sig_pur + sigma_n * randn(1, M);
            
            win = hann(M, 'periodic');
            sig_win = sig_bruite .* win';
            
            P1 = abs(fft(sig_win));
            [valMax, l] = max(P1(1:M/2));
            
            if P1(l+1) > P1(l-1)
                alpha = P1(l+1) / valMax;
                delta_est = (2*alpha - 1) / (alpha + 1);
            else
                alpha = P1(l-1) / valMax;
                delta_est = -(2*alpha - 1) / (alpha + 1);
            end
            deltas_estimes(run) = delta_est;
        end
        var_ipdft_sim(idx_snr) = var(deltas_estimes - delta_k);
    end
    
    % Affichage dans un sous-graphe
    subplot(2, 2, d_idx);
    semilogy(SNR_dB_range, var_crlb_th, 'k--', 'LineWidth', 2); hold on;
    semilogy(SNR_dB_range, var_ipdft_th, 'r-', 'LineWidth', 2);
    semilogy(SNR_dB_range, var_ipdft_sim, 'b.', 'MarkerSize', 12);
    
    grid on;
    title(sprintf('\\delta_k = %.2f', delta_k));
    xlabel('SNR (dB)'); 
    ylabel('Variance (Log)');
    if d_idx == 1
        legend('CRLB (Eq. 14)', 'Théorie IpDFT (Eq. 13)', 'Simulation Pratique');
    end
end
