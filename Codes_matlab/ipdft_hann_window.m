%% INITIALISATION
clear; clc; close all;
fprintf('=== IpDFT AVEC FENÊTRE DE HANNING (H=2) ===\n');

% Paramètres
fs = 1000;              % 1 kHz
N = 1024;               % Taille de la FFT
t = 0:1/fs:(N-1)/fs;    % Axe temps

%% PARTIE 1 : VISUALISATION (Rectangulaire vs Hanning)
% On place une fréquence entre deux bins pour créer de la fuite
f_test = 40.5 * (fs/N); 
sig = cos(2*pi*f_test*t);

% 1. Cas Rectangulaire (Brut)
spec_rect = abs(fft(sig));
spec_rect = spec_rect(1:N/2+1);
spec_rect_dB = 20*log10(spec_rect / max(spec_rect)); % En décibels

% 2. Cas Hanning (Fenêtré)
win = hann(N, 'periodic');  % Création de la fenêtre (H=2)
sig_win = sig .* win';      % Application (Multiplication point par point)
spec_hann = abs(fft(sig_win));
spec_hann = spec_hann(1:N/2+1);
spec_hann_dB = 20*log10(spec_hann / max(spec_hann));

% Affichage comparatif
figure('Name', 'Effet du Fenêtrage (Side Lobe Decay)', 'Color', 'white');
f_axis = fs*(0:(N/2))/N;
plot(f_axis, spec_rect_dB, 'r--', 'LineWidth', 1); hold on;
plot(f_axis, spec_hann_dB, 'b', 'LineWidth', 2);
xlim([0 100]); ylim([-100 0]);
grid on;
title('Comparaison des Spectres : Rectangulaire vs Hanning');
legend('Rectangulaire (Lobes hauts)', 'Hanning (Lobes écrasés - Decay)');
xlabel('Fréquence (Hz)'); ylabel('Amplitude Normalisée (dB)');


%% PARTIE 2 : ESTIMATEUR IpDFT ADAPTÉ (Formule H=2)
fprintf('Lancement du Sweep Test avec formule H=2...\n');

% Plage de test (Autour du Bin 40)
bin_start = 40;
freqs_injectees = (bin_start : 0.05 : bin_start + 2) * (fs/N); 

erreurs_dft = [];
erreurs_ipdft_h2 = [];

for f_vrai = freqs_injectees
    
    % A. Génération Signal
    sig = cos(2*pi*f_vrai*t);
    
    % B. Application Fenêtre Hanning
    % (Indispensable pour utiliser la formule H=2)
    win = hann(N, 'periodic'); 
    sig_win = sig .* win';     
    
    % C. FFT
    Y = fft(sig_win);
    P1 = abs(Y(1:N/2+1));
    
    % D. Recherche du Pic (l)
    [valMax, l] = max(P1);
    
    % E. Calcul du Delta (Formule Spécifique Hanning !)
    delta = 0;
    
    if l > 1 && l < length(P1)
        amp_G = P1(l-1);
        amp_D = P1(l+1);
        
        % La formule change ici par rapport à la fenêtre rectangulaire
        % Formule : delta = (2*alpha - 1) / (alpha + 1)
        
        if amp_D > amp_G
            % Voisin Droite
            alpha = amp_D / valMax;
            delta = (2*alpha - 1) / (alpha + 1);
        else
            % Voisin Gauche
            alpha = amp_G / valMax;
            delta = - (2*alpha - 1) / (alpha + 1); % Signe négatif
        end
    end
    
    % F. Estimation Finale
    idx_corrige = (l - 1) + delta;
    f_estime = idx_corrige * (fs / N);
    
    % Stockage Erreur
    erreurs_ipdft_h2 = [erreurs_ipdft_h2, abs(f_estime - f_vrai)];
    
    % (Comparaison avec DFT classique sur signal non fenêtré pour référence)
    [~, l_raw] = max(abs(fft(sig)));
    f_dft_raw = (l_raw-1)*fs/N;
    erreurs_dft = [erreurs_dft, abs(f_dft_raw - f_vrai)];
end

%% PARTIE 3 : RÉSULTATS DE PERFORMANCE
figure('Name', 'Performance IpDFT Hanning', 'Color', 'white');

subplot(2,1,1);
plot(freqs_injectees, freqs_injectees, 'k--'); hold on;
stairs(freqs_injectees, freqs_injectees + erreurs_dft, 'r');
plot(freqs_injectees, freqs_injectees + erreurs_ipdft_h2, 'b.-');
title('Suivi de Fréquence (Hanning)');
legend('Vraie', 'DFT Standard', 'IpDFT (H=2)'); grid on;

subplot(2,1,2);
plot(freqs_injectees, erreurs_ipdft_h2, 'b', 'LineWidth', 2);
title('Erreur résiduelle avec fenêtre de Hanning');
xlabel('Fréquence (Hz)'); ylabel('Erreur (Hz)');
grid on;

fprintf('Erreur Moyenne IpDFT (Hanning) : %.6f Hz\n', mean(erreurs_ipdft_h2));