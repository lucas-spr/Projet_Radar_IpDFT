%% ÉVALUATION DES PERFORMANCES IPDFT VS CRAMER-RAO (CRLB)

clear; clc; close all;

% Paramètres théoriques
M = 1024;               % Longueur de l'acquisition
fs = 1000;              % Fréquence d'échantillonnage
t = 0:1/fs:(M-1)/fs;    
Ak = 2;                 % Amplitude de la cible (A_k)
bin_k = 13;             % Position (l_k)
delta_k = 0.3;          % Décalage fractionnaire (delta_k)

% Fréquence exacte
fk = (bin_k + delta_k) * (fs / M);

% Plage de SNR testée
SNR_dB_range = 40:10:100;
nb_runs = 1000;         % Nombre de tirages Monte-Carlo par SNR

% Tableaux pour stocker les variances
var_ipdft_simulee = zeros(1, length(SNR_dB_range));
var_crlb_theorique = zeros(1, length(SNR_dB_range));

fprintf('Lancement de la simulation\n');

for idx_snr = 1:length(SNR_dB_range)
    SNR_dB = SNR_dB_range(idx_snr);
    
    % Calcul de la variance du bruit sigma_n^2 à partir du SNR
    % SNR_dB = 10 * log10( (Ak^2 / 2) / sigma_n^2 )
    Puissance_Signal = (Ak^2) / 2;
    SNR_lineaire = 10^(SNR_dB / 10);
    sigma_n_carre = Puissance_Signal / SNR_lineaire;
    
    % Calcul du CRAMER-RAO LOWER BOUND
    % var(delta_CR) = (6 / pi^2) * (sigma_n^2 / (M * Ak^2))
    var_crlb_theorique(idx_snr) = (6 / pi^2) * (sigma_n_carre / (M * Ak^2));
    
    deltas_estimes = zeros(1, nb_runs);
    
    % Boucle de Monte-Carlo (1000 runs)
    for run = 1:nb_runs
        % Signal pur avec phase aléatoire
        phi = rand * 2 * pi;
        sig_pur = Ak * cos(2*pi*fk*t + phi);
        
        % Ajout du bruit blanc gaussien (AWGN)
        bruit = sqrt(sigma_n_carre) * randn(1, M);
        sig_bruite = sig_pur + bruit;
        
        % Application de la fenêtre de Hann (H=2)
        win = hann(M, 'periodic');
        sig_win = sig_bruite .* win';
        
        % FFT et Module
        P1 = abs(fft(sig_win));
        
        % Algorithme IpDFT classique (H=2)
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
    
    % On calcule la variance des erreurs de notre simulation
    erreurs = deltas_estimes - delta_k;
    var_ipdft_simulee(idx_snr) = var(erreurs);
end

%% AFFICHAGE COMPARATIF
figure('Name', 'Comparaison IpDFT vs CRLB', 'Color', 'white');

% On trace en échelle logarithmique (très important pour les variances)
semilogy(SNR_dB_range, var_crlb_theorique, 'k--', 'LineWidth', 2); hold on;
semilogy(SNR_dB_range, var_ipdft_simulee, 'b-o', 'LineWidth', 1.5);

grid on;
title('Variance de l''estimation de \delta_k en fonction du SNR');
xlabel('SNR (dB)');
ylabel('Variance de l''erreur (Echelle Log)');
legend('Borne de Cramer-Rao (Eq. 14)', 'Variance IpDFT (H=2) Simulée');

fprintf('Simulation terminée ! Observez comment la courbe IpDFT suit la borne théorique.\n');