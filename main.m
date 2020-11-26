clc
clear all
warning('off', 'all');

M = 2048;

bb1 = randn(1, M);
bb2 = 9 * randn(1,M) + 7;

%% représentations fréquentielles

% spectrogramme bb1
figure,
subplot(2,1,1);
plot(bb1);
xlim([0;2048]);
subplot(2,1,2);
spectrogram(bb1, 64,0, 'yaxis');
ylabel('Frequences');

% spectrogramme bb2
figure,
subplot(2,1,1);
plot(bb2);
xlim([0;2048]);
subplot(2,1,2);
spectrogram(bb2, 64,0, 'yaxis');
ylabel('Frequences');

%% Etapes DFA

clc
clear all

M = 2001;
N = 64;
L = round(M / N);
bb1 = randn(1, M);
bb2 = 9 * randn(1,M) + 7;

profilY = ProfilY(bb1, M);

% calcul des coefficients qui fit la courbe
coefsX = CalculCoefsProfil(profilY, L, N);

% calcul du profil global
F = CalculProfilGlobal(bb1, L, N);

%% calcul du profil global pour différents N

clc
clear all

M = 2001;
 bb1 = randn(1, M);
bb2 = 9 * randn(1,M) + 7;

Ns = [13, 17, 21, 27, 35, 47, 59, 77, 99, 135, 189];
Fs = zeros(size(Ns));
for i = 1:size(Ns,2)
   N = Ns(i);
   Fs(i) = CalculProfilGlobal(bb1, M, N);
end

% calcul du fit pour obtenir le coefficient de régularité alpha
logF = log(Fs);
logN = log(Ns);
coefsDroite = polyfit(logN, logF, 1);
alpha = coefsDroite(1);
fit = polyval(coefsDroite, logN);
display(alpha);

hold on
plot(logN, logF, 'o');
plot(logN, fit, 'r');
xlabel('log( F(N) )');
ylabel('log( N )');
hold off

%% représentation profil + fit

clc
clear all

M = 2001;
bb1 = randn(1, M);
bb2 = 9 * randn(1,M) + 7;
N = 189;
L = floor(M/N);

profilY = ProfilY(bb1, M);
coefs = CalculCoefsProfil(profilY, L, N);
fit = CalculFit(bb1, M, N, L);


hold on
plot(profilY);
plot(fit, 'r');

for l = 1:L
    xline(l*N);
end
hold off

%% statistiques sur la régularité bb centré

clc
clear all

M = 2001;
alphasCentres = zeros(1,50);
alphasNonCentres = zeros(1,50);
w = waitbar(0, 'traitement des régularités');
for i = 1:50
    bbCentre = randn(1, M);
    bbNonCentre = 9 * randn(1,M) + 7;
    
    alphasCentres(i) = Regularite(bbCentre);
    alphasNonCentres(i) = Regularite(bbNonCentre);
    
    waitbar(i/50);
end
close(w);

mu1 = mean(alphasCentres);
sd1 = std(alphasCentres);
mu2 = mean(alphasNonCentres);
sd2 = std(alphasNonCentres);
%%
hold on 
errorbar([mu1, mu2], [sd1, sd2], 'or');
xlim([0, 3]);
xticks([1,2]);
xticklabels({'Bb centré', 'Bb non centré'});