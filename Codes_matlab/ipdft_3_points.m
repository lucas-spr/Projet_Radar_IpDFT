% SIMULATIONS COMPLÈTES : IpDFT (2 Points) vs WMIpDFT (3 Points)
clear; clc; close all;

% Paramètres globaux du radar et du signal
fs = 1000;              % Fréquence d'échantillonnage
M = 1024;               % Taille de la FFT
t = 0:1/fs:(M-1)/fs;    % Axe temps
res_freq = fs/M;        % Résolution (taille d'un bin)
A = 1;                  % Amplitude

% Fenêtre de Hann (H=2)
win = hann(M, 'periodic');

% SIMULATION 1 : COMPARAISON DU BIAIS SYSTÉMATIQUE (L'Erreur en "S")
fprintf('Lancement Sim 1 : Calcul du Biais Systématique...\n');

bin_k = 15; 
deltas = -0.49:0.01:0.49;

biais_2pts = zeros(1, length(deltas));
biais_3pts = zeros(1, length(deltas));

for i = 1:length(deltas)
    delta_vrai = deltas(i);
    f_vrai = (bin_k + delta_vrai) * res_freq;
    
    % Signal pur
    sig = A * cos(2*pi*f_vrai*t); 
    P1 = abs(fft(sig .* win'));
    
    [valMax, l] = max(P1(1:M/2));
    X_g = P1(l-1); X_c = valMax; X_d = P1(l+1);
    
    %IpDFT 2 Points
    if X_d > X_g
        a2 = X_d / X_c;
        d2 = (2*a2 - 1) / (a2 + 1);
    else
        a2 = X_g / X_c;
        d2 = -(2*a2 - 1) / (a2 + 1);
    end
    biais_2pts(i) = d2 - delta_vrai;
    
    %WMIpDFT 3 Points (H=2, J=1) ---
    a3 = (X_d - X_g) / (X_g + X_d + 2*X_c);
    d3 = 2 * a3;
    biais_3pts(i) = d3 - delta_vrai;
end

figure('Name', 'Sim 1: Biais Systématique', 'Color', 'white');
plot(deltas, biais_2pts, 'r', 'LineWidth', 2); hold on;
plot(deltas, biais_3pts, 'b-', 'LineWidth', 2);
yline(0, 'k--');
grid on;
title('Biais Systématique : 2 Points vs 3 Points');
xlabel('Décalage fractionnaire \delta (Vraie valeur)');
ylabel('Erreur (\delta_{estime} - \delta_{vrai})');
legend('IpDFT Classique (2 points)', 'WMIpDFT (3 points - Écrasé)', 'Location', 'Best');

ylim([-max(abs(biais_2pts))*1.2, max(abs(biais_2pts))*1.2]);

% SIMULATION 2 : COMPARAISON VRAIE FRÉQUENCE vs ESTIMATION (L'Escalier)
fprintf('Lancement Sim 2 : Suivi de fréquence...\n');

bin_start = 40;
freqs_injectees = (bin_start : 0.05 : bin_start + 2) * res_freq; 

err_dft = zeros(1, length(freqs_injectees));
err_2pts = zeros(1, length(freqs_injectees));
err_3pts = zeros(1, length(freqs_injectees));

for i = 1:length(freqs_injectees)
    f_vrai = freqs_injectees(i);
    sig = cos(2*pi*f_vrai*t);
    
    % DFT brute (Rectangulaire)
    [~, l_raw] = max(abs(fft(sig)));
    f_dft = (l_raw-1) * res_freq;
    err_dft(i) = abs(f_dft - f_vrai);
    
    % IpDFT (Fenêtrée)
    P1 = abs(fft(sig .* win'));
    [valMax, l] = max(P1(1:M/2));
    X_g = P1(l-1); X_c = valMax; X_d = P1(l+1);
    
    % 2 Points
    if X_d > X_g
        a2 = X_d / X_c; d2 = (2*a2 - 1) / (a2 + 1);
    else
        a2 = X_g / X_c; d2 = -(2*a2 - 1) / (a2 + 1);
    end
    f_2pts = (l - 1 + d2) * res_freq;
    err_2pts(i) = abs(f_2pts - f_vrai);
    
    % 3 Points
    a3 = (X_d - X_g) / (X_g + X_d + 2*X_c);
    f_3pts = (l - 1 + 2*a3) * res_freq;
    err_3pts(i) = abs(f_3pts - f_vrai);
end

figure('Name', 'Sim 2: Suivi et Erreur', 'Color', 'white');
subplot(2,1,1);
plot(freqs_injectees, freqs_injectees, 'k--'); hold on;
stairs(freqs_injectees, freqs_injectees + err_dft, 'r');
plot(freqs_injectees, freqs_injectees + err_2pts, 'g.-');
plot(freqs_injectees, freqs_injectees + err_3pts, 'b.-');
title('Suivi de Fréquence');
legend('Vérité', 'DFT', 'IpDFT (2 pts)', 'WMIpDFT (3 pts)', 'Location', 'NorthWest'); grid on;


subplot(2,1,2);
semilogy(freqs_injectees, err_dft + eps, 'r-o', 'LineWidth', 1.5); hold on;
semilogy(freqs_injectees, err_2pts + eps, 'g-x', 'LineWidth', 1.5);
semilogy(freqs_injectees, err_3pts + eps, 'b-*', 'LineWidth', 1.5);

title('Erreur Absolue (Échelle Logarithmique)');
xlabel('Fréquence (Hz)'); 
ylabel('Erreur (Hz) - Log');
legend('Erreur DFT (10^{-1})', 'Erreur 2 pts (10^{-5})', 'Erreur 3 pts (10^{-15})', 'Location', 'South');
grid on;

% SIMULATION 3 : ROBUSTESSE À LA PHASE INITIALE
fprintf('Lancement Sim 3 : Robustesse à la phase...\n');

phases = linspace(0, 2*pi, 100); % Balayage de la phase de 0 à 360 degrés
test_deltas = [0.1, 0.3, 0.45];  % Test sur 3 positions différentes dans la case

figure('Name', 'Sim 3: Robustesse Phase', 'Color', 'white');
hold on; grid on;
couleurs = ['r', 'g', 'b'];

for d_idx = 1:length(test_deltas)
    delta_vrai = test_deltas(d_idx);
    f_vrai = (bin_k + delta_vrai) * res_freq;
    erreurs_phase = zeros(1, length(phases));
    
    for p_idx = 1:length(phases)
        phase_actuelle = phases(p_idx);
        
        % Signal avec phase variable
        sig = A * cos(2*pi*f_vrai*t + phase_actuelle);
        P1 = abs(fft(sig .* win'));
        
        [valMax, l] = max(P1(1:M/2));
        X_g = P1(l-1); X_c = valMax; X_d = P1(l+1);
        
        % WMIpDFT (3 Points)
        a3 = (X_d - X_g) / (X_g + X_d + 2*X_c);
        d3 = 2 * a3;
        
        erreurs_phase(p_idx) = d3 - delta_vrai;
    end
    
    plot(phases, erreurs_phase, couleurs(d_idx), 'LineWidth', 2, ...
         'DisplayName', ['\delta = ' num2str(delta_vrai)]);
end

title('Robustesse de la méthode à 3 points (WMIpDFT) face à la phase');
xlabel('Phase initiale (Radians)');
ylabel('Erreur d''estimation (\delta_{est} - \delta_{vrai})');
legend('Location', 'Best');
xlim([0 2*pi]);