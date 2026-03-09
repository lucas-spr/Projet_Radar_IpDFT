%% VISUALISATION DES VALEURS IPDFT (Test de Stabilité)
clear; clc; close all;
fprintf('=== TEST : AFFICHAGE DES VALEURS IPDFT ===\n');

% --- 1. Paramètres ---
fs = 1000;              % Fréquence d'échantillonnage
N = 1024;               % Taille FFT
t = 0:1/fs:(N-1)/fs;    % Axe temps

% On choisit le delta que le système doit trouver
delta_cible = 0.42;     

% On place le signal sur un bin arbitraire (ex: 50)
bin_central = 50;
f_vraie = (bin_central + delta_cible) * (fs/N);

fprintf('Objectif : Le système doit trouver Delta = %.4f\n', delta_cible);
fprintf('Condition : 100 itérations avec PHASE ALÉATOIRE.\n\n');

% Tableaux pour stocker les résultats
nb_tests = 100;
valeurs_ipdft = zeros(1, nb_tests);
erreurs = zeros(1, nb_tests);

%% --- 2. Boucle de Simulation ---
for i = 1:nb_tests
    
    % A. Génération (Phase Aléatoire)
    phi = rand * 2 * pi; 
    sig = cos(2*pi*f_vraie*t + phi);
    
    % B. Fenêtrage Hanning (Important pour la précision)
    win = hann(N, 'periodic');
    sig_win = sig .* win';
    
    % C. FFT
    P1 = abs(fft(sig_win));
    
    % D. IpDFT (Algorithme)
    [valMax, l] = max(P1(1:N/2));
    
    % Logique Voisin (Gauche ou Droite)
    delta_calc = 0;
    
    % Sécurité bords
    if l > 1 && l < length(P1)
        if P1(l+1) > P1(l-1)
            % Voisin Droite
            alpha = P1(l+1) / valMax;
            delta_calc = (2*alpha - 1) / (alpha + 1);
        else
            % Voisin Gauche
            alpha = P1(l-1) / valMax;
            delta_calc = - (2*alpha - 1) / (alpha + 1);
        end
    end
    
    % E. Stockage
    valeurs_ipdft(i) = delta_calc;
    erreurs(i) = delta_calc - delta_cible;
    
    % Affichage textuel pour les 10 premières valeurs seulement
    if i <= 10
        fprintf('Itération %d : Calculé = %.6f (Erreur : %.1e)\n', i, delta_calc, erreurs(i));
    end
end

%% --- 3. Affichage Graphique ---
figure('Name', 'Valeurs IpDFT Brutes', 'Color', 'white');

% GRAPHE 1 : Les valeurs brutes (Ce que vous avez demandé)
subplot(2,1,1);
plot(valeurs_ipdft, 'b.-', 'MarkerSize', 10, 'LineWidth', 0.5); hold on;
yline(delta_cible, 'r--', 'LineWidth', 2); % La ligne cible
title(sprintf('Valeurs de Delta calculées par IpDFT (Cible : %.4f)', delta_cible));
ylabel('Delta Calculé');
legend('Valeur IpDFT (à chaque essai)', 'Cible Théorique');
grid on;

%On zoome fort autour de la valeur
ylim([delta_cible - 0.0001, delta_cible + 0.0001]); 

% GRAPHE 2 : L'Erreur
subplot(2,1,2);
plot(erreurs, 'k', 'LineWidth', 1.5);
title('Erreur de l''estimation (Ecart par rapport à la cible)');
xlabel('Numéro de l''essai (Phase aléatoire)');
ylabel('Erreur');
grid on;
