clc
clear all
warning('off', 'all');

M = 2001;

bb1 = randn(1, M);
bb2 = 9 * randn(1,M) + 7;

%% représentations fréquentielles

close all
% spectrogramme bb1
figure,
subplot(2,1,1);
plot(bb1);
xlabel('Temps en secondes');
ylabel('Intensité signal');
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});

subplot(2,1,2);
spectrogram(bb1, 64,0, 'yaxis');
xlabel('Temps en secondes');
ylabel('Frequences');
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});

% spectrogramme bb2
figure,
subplot(2,1,1);
plot(bb2);
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});

subplot(2,1,2);
spectrogram(bb2, 64,0, 'yaxis');
xlabel('Temps en secondes');
ylabel('Frequences');
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});

%% Etapes DFA

clc
clear all
warning('off', 'all');

M = 2001;
N = 64;
L = round(M / N);
bb1 = randn(1, M);
bb2 = 9 * randn(1,M) + 7;

profilY = ProfilY(bb1, M);

% calcul des coefficients qui fit la courbe
coefsX = CalculCoefsProfil(profilY, L, N);

% calcul du profil global
F = CalculF(bb1, M, N);

for l = 1:L
    for n = 1:N
        m = (l-1) * N + n;
        tempX = polyval(coefsX(l,:), m);
        newY(m, 1) = profilY(m) - tempX;
    end
end

%% calcul du profil global pour différents N

clc
clear all
warning('off', 'all');

M = 2001;
bb1 = randn(1, M);
bb2 = 9 * randn(1,M) + 7;

Ns = [13, 17, 21, 27, 35, 47, 59, 77, 99, 135, 189];
Fs = zeros(size(Ns)); % vecteurs des valeurs F à calculer en fonction de N
for i = 1:size(Ns,2)
   N = Ns(i);
   Fs(i) = CalculF(bb1, M, N);
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
xlabel('log( N )');
ylabel('log( F(N) )');
hold off

%% représentation profil + fit

clc
clear all
warning('off', 'all');

M = 2001;
bb1 = randn(1, M);
bb2 = 9 * randn(1,M) + 7;
N = 189;
L = floor(M/N);

profilY = ProfilY(bb1, M);
coefs = CalculCoefsProfil(profilY, L, N);
fit = CalculFit(bb1, M, N, L);

figure,
subplot(3,1,1), plot(profilY), xlim([0,2001]);
subplot(3,1,2), plot(fit, 'r'), xlim([0,2001]);
subplot(3,1,3),
hold on
plot(profilY), xlim([0,2001]);
plot(fit, 'r');

for l = 1:L
    xline(l*N);
end

hold off

%% DFA statistiques sur la régularité 

clc
clear all
warning('off', 'all');

M = 2001;
alphasCentres = zeros(1,50);
alphasNonCentres = zeros(1,50);
w = waitbar(0, 'traitement des régularités');
for i = 1:50
    bbCentre = randn(1, M);
    bbNonCentre = 9 * randn(1,M) + 7;
    
    alphasCentres(i) = Regularite(bbCentre, 'DFA');
    alphasNonCentres(i) = Regularite(bbNonCentre,'DFA');
    
    waitbar(i/50);
end
close(w);

mu1 = mean(alphasCentres);
sd1 = std(alphasCentres);
mu2 = mean(alphasNonCentres);
sd2 = std(alphasNonCentres);

% figure avec moyennes et écarts types
errorbar([mu1, mu2], [sd1, sd2], 'or');
xlim([0, 3]);
xticks([1,2]);
xticklabels({'Bb centré', 'Bb non centré'});

%% DMA étape par étape

clc
clear all
warning('off', 'all');

M = 2001;
N = 64;
bb1 = randn(1, M);
filtre = ones(1, N) / N;

profil = ProfilY(bb1, M);
tendance = conv(profil, filtre, 'same');

subplot(3,1,1), plot(profil), xlim([0, 2001]);
subplot(3,1,2), plot(tendance, 'r'), xlim([0, 2001]);
subplot(3,1,3), xlim([0, 2001]);
hold on
plot(profil);
plot(tendance, 'r');
hold off

%% représentations du filtre 

clc
clear all
close all
warning('off', 'all');

M = 2001;
N = 64;
bb1 = randn(1, M);
profil = ProfilY(bb1, M);

b = zeros(1, N+1);
b(1) = 1;
b(N+1) = -1;
a = [N, -N];

z = roots(b);
p = roots(a);
figure, zplane(b, a);

[h, w] = freqz(b, a, M, 'whole');
figure, plot(w/pi, 20*log10(abs(h)));
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

[phi, wPhi] = phasez(b, a, M);
figure, plot(wPhi, phi);

%plot de la tendance par rapport au profil
newY = filter(b, a, profil);
retardPhi = (phi(20) - phi(2)) / (wPhi(20) - wPhi(2));
tendance = newY(abs(floor(retardPhi)):size(profil,2));

figure,
hold on
plot(profil);
% plot(newY);
plot(tendance);
xlabel('Temps en secondes');
ylabel('Signal');
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});
hold off

% TF pour voir type de filtre

tf = fft(tendance);
plot(linspace(-0.5, 0.5, size(tf,2)), fftshift(abs(tf)));
xlabel('Fréquences normalisées');
ylabel('Module TF');
% xticks([-0.5,0,0.5]);
%% calcul régularité avec DMA

clc
clear all
warning('off', 'all');
close all

M = 2001;
bb1 = randn(1, M);

[alpha1, beta1, logN1, logF1] = Regularite(bb1, 'DMA');
fit1 = polyval([alpha1, beta1], logN1);
[alpha2, beta2, logN2, logF2] = Regularite(bb1, 'DFA');
fit2 = polyval([alpha2, beta2], logN2);

hold on;
plot(logN1, logF1, 'or');
plot(logN1, fit1);
plot(logN2, logF2, 'ob');
plot(logN2, fit2);
hold off

%% DMA statistiques

clc
clear all
warning('off', 'all');

M = 2001;
alphasCentres = zeros(1,50);
alphasNonCentres = zeros(1,50);
w = waitbar(0, 'traitement des régularités');
for i = 1:50
    bbCentre = randn(1, M);
    bbNonCentre = 9 * randn(1,M) + 7;
    
    alphasCentres(i) = Regularite(bbCentre, 'DMA');
    alphasNonCentres(i) = Regularite(bbNonCentre,'DMA');
    
    waitbar(i/50);
end
close(w);

mu1 = mean(alphasCentres);
sd1 = std(alphasCentres);
mu2 = mean(alphasNonCentres);
sd2 = std(alphasNonCentres);

% figure avec moyennes et écarts types
errorbar([mu1, mu2], [sd1, sd2], 'or');
xlim([0, 3]);
xticks([1,2]);
xticklabels({'Bb centré', 'Bb non centré'});


%% statistiques DFA/DMA

clc
clear all
warning('off', 'all');

M = 2001;
nbAlpha = 100;
alphasDFA = zeros(1,nbAlpha);
alphasDMA = zeros(1,nbAlpha);
w = waitbar(0, 'traitement des régularités');
for i = 1:nbAlpha
    bbCentre = randn(1, M);
    
    alphasDFA(i) = Regularite(bbCentre, 'DFA');
    alphasDMA(i) = Regularite(bbCentre,'DMA');
    
    waitbar(i/nbAlpha);
end
close(w);

mu1 = mean(alphasDFA);
sd1 = std(alphasDFA);
mu2 = mean(alphasDMA);
sd2 = std(alphasDMA);

% figure avec moyennes et écarts types
errorbar([mu1, mu2], [sd1, sd2], 'or');
xlim([0, 3]);
xticks([1,2]);
xticklabels({'DFA', 'DMA'});

%% analyse signaux réels

clc
clear
close all

data = load('dataEEG2020.mat');

%création du dataframe
alpha = zeros(28, 1);
electrode = repmat(' ', 28,1);
method = repmat(' ', 28,1);
df = table(electrode, method, alpha);
% sujet = repmat(' ', 28,1);
% etat = repmat(' ', 28,1);
% df = table(sujet, electrode, etat, alpha);



% analyse avec DFA -------------------------------

nbElectrodes = 2;
nbPhases = size(data.dataEEG2020e7, 1);
nbSujets = size(data.dataEEG2020e7, 2);

% DFA electrode 1
alphaDFAe1 = zeros(1, nbPhases*nbSujets);
index = 1;
for p = 1:nbPhases
    for s = 1:nbSujets
        
        y = cell2mat(data.dataEEG2020e7(p,s))';
        alphaDFAe1(1, index) = Regularite(y, 'DFA');
        index = index +1;

    end
    
end

% DFA electrode 2
alphaDFAe2 = zeros(1, nbPhases*nbSujets);
index = 1;
for p = 1:nbPhases
    for s = 1:nbSujets
        
        % le signal est à transposer car on l'obtient en vecteur colonne
        y = cell2mat(data.dataEEG2020e8(p,s))'; 
        alphaDFAe2(1, index) = Regularite(y, 'DFA');
        index = index +1;

    end
    
end



% analyse avec DMA -------------------------------

nbElectrodes = 2;
nbPhases = size(data.dataEEG2020e7, 1);
nbSujets = size(data.dataEEG2020e7, 2);

% DFA electrode 1
alphaDMAe1 = zeros(1, nbPhases*nbSujets);
index = 1;
for p = 1:nbPhases
    for s = 1:nbSujets
        
        y = cell2mat(data.dataEEG2020e7(p,s))';
        alphaDMAe1(1, index) = Regularite(y, 'DMA');
        index = index +1;
        
    end
    
end

% DFA electrode 2
alphaDMAe2 = zeros(1, nbPhases*nbSujets);
index = 1;
for p = 1:nbPhases
    for s = 1:nbSujets
        
        % le signal est à transposer car on l'obtient en vecteur colonne
        y = cell2mat(data.dataEEG2020e8(p,s))'; 
        alphaDMAe2(1, index) = Regularite(y, 'DMA');
        index = index +1;
        
    end
    
end

% Figures -----------------------------------------------

% plot DFA
figure,
hold on
plot(alphaDFAe1);
plot(alphaDFAe2, 'r');
legend('e1', 'e2');
title('alpha obtenus par DFA');
hold off

% plot DMA
figure,
hold on
plot(alphaDMAe1);
plot(alphaDMAe2, 'r');
legend('e1', 'e2');
title('alpha obtenus par DMA');
hold off

% plot des stats 
moyennes = [mean(alphaDFAe1), mean(alphaDFAe2), mean(alphaDMAe1), mean(alphaDMAe2)];
sd = [std(alphaDFAe1), std(alphaDFAe2), std(alphaDMAe1), std(alphaDMAe2)];

figure,
errorbar(moyennes, sd, 'or');
xlim([0,5]);
xticks([1,2,3,4]);
xticklabels({'DFAe1', 'DFAe2','DMAe1','DMAe2'});










