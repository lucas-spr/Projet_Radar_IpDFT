%Déclaration des paramètres
% Déclaration des paramètres
A0 = 0; % Amplitude de la composante continue
A1 = 1; % Amplitude de la composante sinusoïdale
fech = 1000; % Fréquence d'échantillonnage en Hz
M = 2048; % Nombre de points pour la FFT
N=1000;
m=0:N-1;
k0=50;
f1=k0*fech/N;

%Création d'un signal x
x = sin(2*pi*f1*m/fech);


%Affichage de sa transformée de Fourier décalée
X = fft(x, M); % Calcul de la FFT et décalage
Xshift=fftshift(X);

f = (-M/2:M/2-1)*(fech/M); % Fréquence correspondante
figure;
plot(f, abs(Xshift)); % Affichage du module de la FFT
xlabel('Fréquence (Hz)');
ylabel('Amplitude');
title('Transformée de Fourier du signal x');


%Calcul du maximum de la transformée de Fourier
[maxAmplitude, maxIndex] = max(abs(Xshift)); % Calcul du maximum et de son index
maxFrequency = f(maxIndex); % Fréquence correspondante au maximum