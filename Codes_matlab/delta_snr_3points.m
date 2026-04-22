%% ÉVALUATION STATISTIQUE WMIpDFT (3 POINTS) SOUS BRUIT AWGN
clear; clc; close all;

fprintf('=== SIMULATION DE MONTE-CARLO : WMIpDFT (3 POINTS) ===\n');

% Paramètres de la simulation
M = 1024;                   % Taille de la FFT
fs = 1000;                  % Fréquence d'échantillonnage
t = 0:1/fs:(M-1)/fs;        % Axe temps
A = 1;                      % Amplitude du signal
bin_k = 15;                 % Bin de référence
win = hann(M, 'periodic');  % Fenêtre de Hanning (H=2)

% Paramètres de l'étude
deltas_test = [0.10, -0.40, 0.30, -0.25];
snr_db = 40:10:100;         % SNR de 40 dB à 100 dB
N_runs = 1000;              % Nombre de tirages Monte-Carlo par point

% Préparation de la figure (4 subplots)
figure('Name', 'Variance WMIpDFT 3 points', 'Color', 'white', 'Position', [100, 100, 1000, 800]);

for d_idx = 1:length(deltas_test)
    delta_vrai = deltas_test(d_idx);
    f_vrai = (bin_k + delta_vrai) * (fs / M);
    
    var_sim = zeros(1, length(snr_db));
    crlb = zeros(1, length(snr_db));
    
    fprintf('Calcul pour delta = %.2f...\n', delta_vrai);
    
    for s_idx = 1:length(snr_db)
        snr_actuel = snr_db(s_idx);
        
        % Calcul de l'écart-type du bruit temporel à partir du SNR
        % SNR_dB = 10*log10( P_signal / P_bruit )
        % P_signal pour un cosinus = A^2 / 2
        P_sig = (A^2) / 2;
        P_noise = P_sig / (10^(snr_actuel / 10));
        sigma_n = sqrt(P_noise);
        
        % Calcul de la Borne de Cramer-Rao (Eq. 14)
        crlb(s_idx) = (6 / pi^2) * (sigma_n^2 / (M * A^2));
        
        erreurs_run = zeros(1, N_runs);
        
        %MONTE-CARLO
        for run = 1:N_runs
            phase = rand() * 2 * pi;
            
            sig_pur = A * cos(2*pi*f_vrai*t + phase);
            bruit = sigma_n * randn(1, M);
            sig_bruite = sig_pur + bruit;
            
            % FFT
            P1 = abs(fft(sig_bruite .* win'));
            

            l = bin_k + 1; 
            
            % Extraction des 3 points (Voisin Gauche, Centre, Voisin Droit)
            X_g = P1(l-1);
            X_c = P1(l);
            X_d = P1(l+1);
            
            % Application de la formule WMIpDFT (3 points, H=2)
            alpha_3 = (X_d - X_g) / (X_g + X_d + 2*X_c);
            delta_estime = 2 * alpha_3;
            
            erreurs_run(run) = delta_estime - delta_vrai;
        end
        
        % Calcul de la variance expérimentale sur les 1000 tirages
        var_sim(s_idx) = var(erreurs_run);
    end
    
    %AFFICHAGE DU SUBPLOT
    subplot(2, 2, d_idx);
    semilogy(snr_db, crlb, 'k--', 'LineWidth', 2); hold on;
    semilogy(snr_db, var_sim, 'b.', 'MarkerSize', 15);
    
    title(sprintf('\\delta_k = %.2f', delta_vrai));
    xlabel('SNR (dB)');
    ylabel('Variance (Log)');
    grid on;
    
    ylim([10^-16, 10^-6]); 
    xlim([40, 100]);
    
    if d_idx == 1
        legend('CRLB (Eq. 14)', 'Simulation Pratique 3 pts', 'Location', 'NorthEast');
    end
end

fprintf('Terminé !\n');