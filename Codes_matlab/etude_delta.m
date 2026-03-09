%% ÉTUDE DE L'IMPACT DU DELTA (DÉCALAGE)
clear; clc; close all;

% Paramètres de simulation (Petits pour bien voir les barres)
fs = 100;               % Hz
N = 32;                 % 32 points seulement (gros pixels)
t = 0:1/fs:(N-1)/fs;    % Axe temps

% Bin central pour le test (ex: Bin #5)
bin_cible = 5;
res_freq = fs/N;        % Résolution (3.125 Hz)

% Liste des Deltas à tester :
% 0    = Pile sur le bin (Aucune fuite)
% 0.25 = Un petit décalage (Fuite légère)
% 0.5  = Pire cas (Exactement entre deux bins)
deltas_test = [0, 0.25, 0.5];

figure('Name', 'Comparaison des Deltas', 'Color', 'white');

% Boucle pour générer et tracer les 3 cas
for i = 1:length(deltas_test)
    d = deltas_test(i);
    
    % 1. Calcul de la fréquence exacte correspondant à ce delta
    % Fréquence = (Index + Delta) * Résolution
    f_test = (bin_cible + d) * res_freq;
    
    % 2. Génération du signal (Fenêtre Rectangulaire pour bien voir la fuite)
    sig = cos(2*pi*f_test*t);
    
    % 3. FFT
    Spec = abs(fft(sig)/N);
    Spec = Spec(1:N/2+1); % On garde la partie positive
    
    % 4. Estimation IpDFT (Calcul Inverse pour vérifier)
    [valMax, l] = max(Spec);
    % Voisin de droite (car nos deltas sont positifs)
    if l < length(Spec)
        alpha = Spec(l+1) / valMax;
        delta_estime = alpha / (1 + alpha); % Formule Rectangulaire
    else
        delta_estime = 0;
    end
    
    % 5. Affichage (Subplot)
    subplot(3, 1, i);
    stem(0:length(Spec)-1, Spec, 'filled', 'LineWidth', 1.5); hold on;
    
    % Ligne verte pour montrer la VRAIE position
    xline(bin_cible + d, 'g--', 'LineWidth', 2);
    
    % Mise en forme
    grid on; xlim([0 10]); ylim([0 0.6]);
    title(sprintf('Cas Delta = %.2f (Freq = %.2f Hz)', d, f_test));
    xlabel('Index du Bin'); ylabel('Amplitude');
    
    % Affichage console
    fprintf('--- Cas Delta = %.2f ---\n', d);
    fprintf('   Pic principal détecté au bin : %d\n', l-1);
    fprintf('   Amplitude Pic (Bin %d)       : %.4f\n', l-1, Spec(l));
    fprintf('   Amplitude Voisin (Bin %d)    : %.4f\n', l, Spec(l+1));
    fprintf('   Delta Estimé par IpDFT       : %.4f (Erreur: %.1e)\n\n', delta_estime, abs(d - delta_estime));
end