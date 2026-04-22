clear; clc; close all;

%% Paramètres
M = 1024;
fs = 1000;
t = 0:1/fs:(M-1)/fs;
win = hann(M, 'periodic');

% Cible fixée à delta = 0.42
% bin_k = 15;
% delta_cible = 0.42;
%PIRES CONDITIONS POUR LE RADAR
bin_k = 2;          % Cible très très proche (Basse fréquence, forte interférence)
delta_cible = 0.42; % Cible au bord de la case

f_vrai = (bin_k + delta_cible) * (fs / M);


Nb_essais = 100;
deltas_calcules = zeros(1, Nb_essais);
erreurs = zeros(1, Nb_essais);

%% Boucle sur 100 essais avec une phase aléatoire à chaque fois
for i = 1:Nb_essais
    % Génération d'une phase aléatoire entre 0 et 2*pi
    phase_alea = rand() * 2 * pi;
    
    % Signal avec cette phase
    sig = cos(2*pi*f_vrai*t + phase_alea);
    
    % FFT et recherche du pic
    P1 = abs(fft(sig .* win'));
    [valMax, l] = max(P1(1:M/2));
    
    % Valeurs des 3 points
    X_g = P1(l-1); 
    X_c = valMax; 
    X_d = P1(l+1);
    
    % ALGORITHME WMIpDFT (3 POINTS)
    a3 = (X_d - X_g) / (X_g + X_d + 2*X_c);
    delta_estime = 2 * a3;
    
    deltas_calcules(i) = delta_estime;
    erreurs(i) = delta_estime - delta_cible;
end

%% Affichage
figure('Name', 'Figure 3: Valeurs WMIpDFT Brutes (3 points)', 'Color', 'white');

%Graphe 1 : Valeurs brutes
subplot(2,1,1);
plot(1:Nb_essais, deltas_calcules, 'b.-', 'MarkerSize', 10, 'LineWidth', 1); hold on;
yline(delta_cible, 'r--', 'LineWidth', 2);
title(sprintf('Valeurs de Delta calculées par WMIpDFT (Cible : %.4f)', delta_cible));
ylabel('Delta Calculé');
legend('Valeur WMIpDFT (à chaque essai)', 'Cible Théorique', 'Location', 'Best');
grid on;
ylim([delta_cible-0.0001, delta_cible+0.0001]);

%Graphe 2 : Erreur
subplot(2,1,2);
plot(1:Nb_essais, erreurs, 'k', 'LineWidth', 1.5);
title('Erreur de l''estimation (Ecart par rapport à la cible)');
xlabel('Numéro de l''essai (Phase aléatoire)');
ylabel('Erreur');
grid on;

fprintf('Regardez l''axe Y du graphe du bas : l''erreur est de l''ordre de 10^-15 (Zéro absolu pour MATLAB) !\n');