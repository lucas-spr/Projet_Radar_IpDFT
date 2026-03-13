clear; clc; close all;
fprintf('=== SIMULATION IPDFT : COMPARAISON THÉORIE / PRATIQUE ===\n');

%% 1. Paramètres
M = 1024;               % Longueur de l'acquisition (Nombre de points)
fs = 1000;              % Fréquence d'échantillonnage (Hz)
t = 0:1/fs:(M-1)/fs;    % Vecteur temps

% Paramètres du signal (Cible unique pour isoler la variance du bruit)
Ak = 2;                 % Amplitude de la cible
bin_k = 13;             % Position entière (l_k)
delta_k = 0.3;          % Décalage fractionnaire (delta_k)
fk = (bin_k + delta_k) * (fs / M); % Fréquence réelle

H = 2;                  % Fenêtre de Hann (Maximum Sidelobe Decay d'ordre 2)

% Paramètres de la simulation
SNR_dB_range = 40:10:100; % Plage de SNR (40 dB à 100 dB)
nb_runs = 1000;           % Itérations Monte-Carlo par point

% Initialisation des vecteurs de résultats
var_crlb_theorique = zeros(1, length(SNR_dB_range));
var_ipdft_theorique = zeros(1, length(SNR_dB_range));
var_ipdft_simulee = zeros(1, length(SNR_dB_range));

fprintf('Lancement de la simulation sur %d itérations par SNR...\n', nb_runs);

%% 2. Boucle Principale
for idx_snr = 1:length(SNR_dB_range)
    SNR_dB = SNR_dB_range(idx_snr);
    
    % Calcul de la variance du bruit (sigma_n^2) à partir du SNR
    Puissance_Signal = (Ak^2) / 2;
    SNR_lineaire = 10^(SNR_dB / 10);
    sigma_n_carre = Puissance_Signal / SNR_lineaire;
    sigma_n = sqrt(sigma_n_carre);
    
    % --- Borne de Cramer-Rao (CRLB) - Équation 14 ---
    var_crlb_theorique(idx_snr) = (6 / pi^2) * (sigma_n_carre / (M * Ak^2));
    
    % --- Variance théorique de l'IpDFT - Équation 13 ---
    var_ipdft_theorique(idx_snr) = calculer_variance_eq13(delta_k, Ak, M, sigma_n, H);
    
    % --- Simulation de l'estimateur IpDFT ---
    deltas_estimes = zeros(1, nb_runs);
    
    for run = 1:nb_runs
        % Signal avec phase aléatoire et ajout de bruit blanc gaussien (AWGN)
        phi = rand * 2 * pi;
        sig_pur = Ak * cos(2*pi*fk*t + phi);
        bruit = sigma_n * randn(1, M);
        sig_bruite = sig_pur + bruit;
        
        % Application de la fenêtre de Hann
        win = hann(M, 'periodic');
        sig_win = sig_bruite .* win';
        
        % FFT et Algorithme IpDFT (Formule H=2)
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
    
    % Variance empirique obtenue sur les 1000 tirages
    erreurs = deltas_estimes - delta_k;
    var_ipdft_simulee(idx_snr) = var(erreurs);
end

%% 3. Affichage Graphique
figure('Name', 'Validation Théorique IpDFT', 'Color', 'white');

% On trace en échelle logarithmique (essentiel pour les variances)
semilogy(SNR_dB_range, var_crlb_theorique, 'k--', 'LineWidth', 2); hold on;
semilogy(SNR_dB_range, var_ipdft_theorique, 'r-', 'LineWidth', 2);
semilogy(SNR_dB_range, var_ipdft_simulee, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);

grid on;
title('Performance de l''IpDFT : Pratique vs Théorie (Éq. 13 et 14)');
xlabel('Rapport Signal sur Bruit (SNR) en dB');
ylabel('Variance de l''erreur')