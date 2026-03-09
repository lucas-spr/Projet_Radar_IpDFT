clear; clc; close all;
fprintf('=== RADAR FMCW : ESTIMATION DU BATTEMENT PAR IPDFT ===\n\n');

%% 1. Paramètres du Radar (Chirp)
fs = 1000;              % Fréquence d'échantillonnage (1 kHz)
T  = 1.0;               % Durée du tir (1 seconde)
B  = 400;               % Bande passante (Sweep de 0 à 400 Hz)
S  = B / T;             % Pente (Slope) du Chirp (400 Hz/s)

t  = 0:1/fs:T-1/fs;     % Axe temps
N  = length(t);         % Nombre de points (1000)
res_freq = fs/N;        % Résolution DFT = 1 Hz

%% 2. La Cible 
% On choisit un retard qui va donner une fréquence fractionnaire
tau_vrai = 0.1042;      % Retard de 104.2 ms
fb_theorique = S * tau_vrai; % Fréquence de battement = 41.68 Hz

fprintf('--- DONNÉES THÉORIQUES ---\n');
fprintf('Retard cible     : %.4f s\n', tau_vrai);
fprintf('Fréquence Batt.  : %.4f Hz\n\n', fb_theorique);

%% 3. Génération des Signaux physiques
% Émission (Tx) : Phase = pi * S * t^2
Tx = cos(pi * S * t.^2);

% Réception (Rx) : Signal retardé
t_rx = t - tau_vrai;
Rx = cos(pi * S * t_rx.^2);
Rx(t < tau_vrai) = 0; % Pas de signal avant l'arrivée de l'écho

% Mélangeur (Mixer) : Multiplication Tx * Rx
Mix = Tx .* Rx;

%% 4. Analyse par IpDFT (Fenêtre de Hanning H=2)
fprintf('--- TRAITEMENT IPDFT ---\n');

% A. Fenêtrage (Crucial pour écraser les lobes des hautes fréquences)
win = hann(N, 'periodic');
Mix_win = Mix .* win';

% B. Calcul de la FFT
P1 = abs(fft(Mix_win));

% C. Recherche du pic basse fréquence (On cherche dans la 1ère moitié)
[valMax, l] = max(P1(1:floor(N/2)));

% D. Interpolation (IpDFT avec formule Hanning)
if P1(l+1) > P1(l-1)
    alpha = P1(l+1) / valMax;
    delta_calc = (2*alpha - 1) / (alpha + 1);
else
    alpha = P1(l-1) / valMax;
    delta_calc = - (2*alpha - 1) / (alpha + 1);
end

% E. Calcul final
idx_estime = (l - 1) + delta_calc;
fb_estime = idx_estime * res_freq;

% F. Déduction du retard (L'objectif final du Radar !)
tau_estime = fb_estime / S;

fprintf('Fréquence DFT brute : %d.0000 Hz (Erreur : %.2f Hz)\n', l-1, abs((l-1)-fb_theorique));
fprintf('Fréquence IpDFT     : %.4f Hz (Erreur : %.4f Hz)\n', fb_estime, abs(fb_estime-fb_theorique));
fprintf('Retard Radar estimé : %.4f s\n\n', tau_estime);

%% 5. Visualisation
figure('Name', 'Traitement du signal Radar', 'Color', 'white');

% Graphe 1 : Le signal de battement dans le temps
subplot(2,1,1);
plot(t, Mix, 'k');
xlim([tau_vrai, tau_vrai+0.2]); % On zoome un peu après le retard
title('Signal en sortie du mélangeur (Battement basse fréquence + bruit HF)');
xlabel('Temps (s)'); ylabel('Amplitude'); grid on;

% Graphe 2 : Le spectre et l'estimation
subplot(2,1,2);
f_axis = fs*(0:(floor(N/2)-1))/N;
% On trace le spectre en zoomant sur les basses fréquences
plot(f_axis, P1(1:floor(N/2)), 'b', 'LineWidth', 1.5); hold on;

% Marqueurs pour montrer l'action de l'IpDFT
xline(l-1, 'r--', 'DFT Standard'); 
xline(fb_estime, 'g-', 'Estimation IpDFT', 'LineWidth', 2);
xline(fb_theorique, 'k:', 'Vérité');

xlim([30 55]); % Zoom autour du pic 41.68 Hz
title('Spectre du battement et localisation fine de la cible');
xlabel('Fréquence (Hz)'); ylabel('Amplitude');
legend('Spectre du battement', 'Case DFT (Erreur)', 'IpDFT (Précis)', 'Valeur Théorique');
grid on;