%% 1. Paramètres "Visuels" (Échelle Audio)
fs = 1000;              % Fréquence d'échantillonnage (1 kHz)
T = 1;                  % Durée du tir (1 seconde)
t = 0:1/fs:T-1/fs;      % Axe temps (1000 points)

B = 400;                % Bande passante (Sweep 0 -> 400 Hz)
S = B / T;              % Pente (Slope) = 400 Hz/s

% Paramètre de la cible (On définit directement le retard tau)
tau = 0.1;              % Retard de l'écho (0.1s)
% Note: Avec la vitesse lumière, cela correspondrait à une cible à 15 000 km !
% Mais pour la simulation mathématique, c'est parfait pour voir les effets.

%% 2. Génération des Signaux
% Signal Émis (Tx)
Tx = cos(pi * S * t.^2);

% Signal Reçu (Rx) - Retardé
% On crée le signal retardé en décalant l'axe temps
t_lag = t - tau;
Rx = cos(pi * S * t_lag.^2);
% Astuce simulation: On met à 0 la partie où l'écho n'est pas encore revenu
Rx(t < tau) = 0; 

%% 3. Le Mélange (Signal de Battement)
% Mathématiquement: Tx * Rx
Mix = Tx .* Rx;

% Le mélange crée deux fréquences : (Tx-Rx) et (Tx+Rx).
% On veut garder uniquement la différence (Basse Fréquence).
% Un simple filtre passe-bas ou l'étude spectrale le montrera.
% La fréquence théorique du battement est :
f_beat_theorique = S * tau; % Doit donner 40 Hz (400 * 0.1)

%% 4. Visualisation
figure('Name', 'Simulation RADAR - Affichage Corrigé', 'Color', 'white');

% A. Signaux Temporels
subplot(3,1,1);
% On convertit le temps en microsecondes pour l'affichage si c'est très court
if T < 1e-3
    scale_t = 1e6;
    unit_t = '(\mus)';
else
    scale_t = 1;
    unit_t = '(s)';
end

plot(t*scale_t, Tx, 'b', t*scale_t, Rx, 'r--');
% COMMANDE IMPORTANTE : On laisse Matlab ajuster le zoom tout seul
axis tight; 
title('Signaux Tx (Bleu) et Rx (Rouge)');
xlabel(['Temps ' unit_t]); grid on;
legend('Emis', 'Reçu (Retardé)');

% B. Spectrogramme (Pour voir les Chirps)
subplot(3,1,2);
% Fenêtre de 128 points (suffisant car signal lent), overlap important
spectrogram(Tx + Rx, 128, 120, 128, fs, 'yaxis'); 
title('Spectrogramme combiné : On voit les deux lignes parallèles !');

% C. Spectre du Signal de Battement (Ce que voit l''IpDFT)
subplot(3,1,3);
L = length(Mix);
P2 = abs(fft(Mix)/L);
P1 = P2(1:L/2+1);
f_axis = fs*(0:(L/2))/L;

plot(f_axis, P1, 'k', 'LineWidth', 1.5);
xlim([0 100]); % On regarde les basses fréquences
title(['Spectre du Battement (Pic attendu vers ' num2str(f_beat_theorique) ' Hz)']);
xlabel('Fréquence (Hz)'); grid on;