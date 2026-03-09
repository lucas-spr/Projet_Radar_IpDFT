 %% ÉTAPE 0 : Illustration du Problème (La Fuite Spectrale)

% Paramètres simplifiés
fs_demo = 100;          % Freq. échantillonnage faible
N_demo = 32;            % Peu de points pour bien voir les "casiers" (bins)
t_demo = 0:1/fs_demo:(N_demo-1)/fs_demo;

% Résolution de la DFT (Largeur d'un bin)
delta_f = fs_demo / N_demo; 
fprintf('Résolution de la DFT (Taille d''un bin) : %.2f Hz\n', delta_f);

% --- CAS A : Le monde idéal (Coherent Sampling) ---
% On choisit une fréquence qui tombe PILE sur le Bin n°4
f_ideal = 4 * delta_f; 
sig_ideal = cos(2*pi*f_ideal*t_demo);

% --- CAS B : La réalité (Spectral Leakage) ---
% On choisit une fréquence qui tombe ENTRE le Bin n°4 et n°5
f_reel = 4.5 * delta_f; 
sig_reel = cos(2*pi*f_reel*t_demo);

% Calcul des FFT
Spec_A = abs(fft(sig_ideal)/N_demo); % Normalisé
Spec_B = abs(fft(sig_reel)/N_demo);

% --- Visualisation Comparative ---
figure('Name', 'Pourquoi IpDFT ?', 'Color', 'white');

% Plot 1 : Cas Idéal
subplot(2,1,1);
stem(0:N_demo-1, Spec_A, 'b', 'LineWidth', 2, 'MarkerFaceColor', 'b');
title(['CAS A : La fréquence (' num2str(f_ideal) ' Hz) est un multiple entier. Pic net.']);
xlabel('Index du Bin'); ylabel('Amplitude'); grid on;
xlim([0 10]); % Zoom sur le début

% Plot 2 : Cas Réel (Fuite)
subplot(2,1,2);
stem(0:N_demo-1, Spec_B, 'r', 'LineWidth', 2, 'MarkerFaceColor', 'r');
hold on; 
% On dessine une ligne verte là où est la VRAIE fréquence
xline(4.5, 'g--', 'Vraie Fréquence', 'LineWidth', 1.5);
title(['CAS B : La fréquence (' num2str(f_reel) ' Hz) est entre deux bins. L''énergie bave (Fuite) !']);
xlabel('Index du Bin'); ylabel('Amplitude'); grid on;
xlim([0 10]);

fprintf('Appuyez sur une touche pour lancer la suite de la simulation...\n');
pause;
close all; 




%% Initialisation
fs = 1000;              % 1 kHz
N = 1024;               % Taille du signal
t = 0:1/fs:(N-1)/fs;    % Vecteur temps

% On va tester des fréquences entre le Bin 40 et le Bin 41
bin_start = 40;
freqs_test = (bin_start : 0.05 : bin_start + 2) * (fs/N); 
% On teste sur 2 bins complets avec des pas très fins

% Tableaux pour stocker les erreurs
erreurs_dft = [];
erreurs_ipdft = [];

fprintf('Début du test de validation...\n');
%Affichage d'une transformée de Fourier qui tombe pile au milieu d'un bin
%pour illustrer la fuite spectralre
signal_test=cos(2*pi*1.5*t);

%% Boucle de Test (Monte Carlo déterministe)
for f_vrai = freqs_test
    
    % 1. Génération du signal (Pur, sans bruit pour valider la formule)
    sig = cos(2*pi*f_vrai*t);
    
    % 2. Calcul FFT
    Y = fft(sig);
    P1 = abs(Y(1:N/2+1));
    
    % 3. Estimateur DFT Classique (Max)
    [valMax, indexMax] = max(P1);
    l = indexMax; 
    % Rappel: index 1 = 0Hz. Freq = (index-1)*fs/N
    f_dft = (l - 1) * (fs / N);
    
    % 4. Estimateur IpDFT
    % Ici on implémente la version robuste, on regarde le pic de gauche et
    % celui de droite
    if l > 1 && l < length(P1)
        amp_G = P1(l-1);
        amp_D = P1(l+1);
        if amp_D > amp_G
            alpha = amp_D / valMax;
            signe = 1;
        else
            alpha = amp_G / valMax;
            signe = -1;
        end
        delta = signe * alpha / (1 + alpha);
    else
        delta = 0; % Sécurité bords
    end
    
    idx_interp = (l - 1) + delta;
    f_ipdft = idx_interp * (fs / N);
    % --- Fin Algo IpDFT ---
    
    % 5. Stockage des erreurs (Valeur Absolue)
    erreurs_dft = [erreurs_dft, abs(f_dft - f_vrai)];
    erreurs_ipdft = [erreurs_ipdft, abs(f_ipdft - f_vrai)];
end

%% Visualisation de la Preuve
figure('Name', 'Preuve de fonctionnement IpDFT', 'Color', 'white');

% Graphe 1 : Comparaison des estimations
subplot(2,1,1);
plot(freqs_test, freqs_test, 'k--', 'LineWidth', 1); hold on;
stairs(freqs_test, freqs_test + erreurs_dft, 'r', 'LineWidth', 1.5); % DFT en escalier
plot(freqs_test, freqs_test + erreurs_ipdft, 'b.-', 'MarkerSize', 8);
title('Comparaison : Vraie Fréquence vs Estimations');
legend('Vraie Fréquence', 'DFT Classique (Escalier)', 'IpDFT (Interpolée)');
grid on;

% Graphe 2 : L'erreur pure
subplot(2,1,2);
plot(freqs_test, erreurs_dft, 'r-o', 'MarkerSize', 4); hold on;
plot(freqs_test, erreurs_ipdft, 'b-o', 'MarkerSize', 4, 'LineWidth', 1.5);
title('Erreur absolue de l''estimation');
xlabel('Fréquence injectée (Hz)');
ylabel('Erreur (Hz)');
legend('Erreur DFT', 'Erreur IpDFT (quasi nulle)');
grid on;

% Affichage quantitatif
rmse_dft = sqrt(mean(erreurs_dft.^2));
rmse_ipdft = sqrt(mean(erreurs_ipdft.^2));
fprintf('\n--- RÉSULTATS ---\n');
fprintf('Erreur Moyenne DFT   : %.5f Hz\n', mean(erreurs_dft));
fprintf('Erreur Moyenne IpDFT : %.5f Hz\n', mean(erreurs_ipdft));
fprintf('GAIN DE PRÉCISION    : x%.1f\n', mean(erreurs_dft)/mean(erreurs_ipdft));