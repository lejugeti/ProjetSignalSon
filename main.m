clc
clear all
warning('off', 'all');

M = 2001;

bb1 = randn(1, M);
bb2 = 9 * randn(1,M) + 7;

%% plot bb

plot(bb1);
xlabel('Temps (s)');
ylabel('A');
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});
title("Représentation bruit blanc centré");
%% tf bb

tf = fft(bb1(:));

plot(fftshift(abs(tf)));
xlim([0, 2001]);
xticks([0, 1000, 2001]);
xticklabels(["-0.5", "0", "0.5"]);
ylabel("|TF(y)|");
xlabel("Fréquences (Hz)");

%% Spectrogrammes

clc
clear all
warning('off', 'all');
close all

M = 2001;

bb1 = randn(1, M);
bb2 = 9 * randn(1,M) + 7;

% spectrogramme bb1
figure, 
spectrogram(bb1, sqrt(M), 'yaxis');
xlabel('Temps en secondes');
ylabel('Frequences');
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});
title("Bruit blanc centré");


% spectrogramme bb2
figure, 
spectrogram(bb2 - mean(bb2), sqrt(M), 'yaxis');
xlabel('Temps en secondes');
ylabel('Frequences');
xlim([0;M]);
xticks([0,500,1000,1500,2000]);
xlabel('Temps en secondes');
xticklabels({'0', '0.5', '1', '1.5', '2'});
title("Bruit blanc non-centré");

% spectrogramme bb2 recentré
figure, 
spectrogram(bb2 - mean(bb2), sqrt(M), 'yaxis');
xlabel('Temps en secondes');
ylabel('Frequences');
xlim([0;M]);
xticks([0,500,1000,1500,2000]);
xlabel('Temps en secondes');
xticklabels({'0', '0.5', '1', '1.5', '2'});
title("Bruit blanc recentré");
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
title("Représentation des log de F par DFA");
text(logN(9), logF(9), sprintf("   alpha = %.3f", alpha));
hold off

%% représentation profil + fit

clc
clear all
warning('off', 'all');

M = 2001;
bb1 = randn(1, M);
bb2 = 9 * randn(1,M) + 7;
N = 64;
L = floor(M/N);

profilY = ProfilY(bb1, M);
coefs = CalculCoefsProfil(profilY, L, N);
fit = CalculFit(bb1, M, N, L);

figure,

%profil
subplot(3,1,1), plot(profilY), xlim([0,2001]);
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});
ylabel("A");
xlabel("Temps (s)");

% fit
subplot(3,1,2), plot(fit, 'r'), xlim([0,2001]);
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});
ylabel("A");
xlabel("Temps (s)");

% profil + fit
subplot(3,1,3),
hold on
plot(profilY), xlim([0,2001]);
plot(fit, 'r');

for l = 1:L
    xline(l*N);
end

xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});
ylabel("A");
xlabel("Temps (s)");
hold off

% comparaison fit pour N différents

profil = ProfilY(bb1, M);
fit1 = CalculFit(bb1, M, 65, floor(M/65));
fit2 = CalculFit(bb1, M, 189, floor(M/189));

figure,
subplot(2,1,1),
hold on
plot(profil), xlim([0,2001]);
plot(fit1, 'r');

xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});
ylabel("A");
xlabel("Temps (s)");
title("Tendance avec N = 65");
hold off

subplot(2,1,2),
hold on
plot(profil), xlim([0,2001]);
plot(fit2, 'r');

xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});
ylabel("A");
xlabel("Temps (s)");
title("Tendance avec N = 189");
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
ylabel("Alpha");
title("Alphas moyens par DFA");
%% DMA étape par étape

clc
clear all
warning('off', 'all');
close all

M = 2001;
N = 64;
bb1 = randn(1, M);

b = zeros(1, N+1);
b(1) = 1;
b(N+1) = -1;
a = [N, -N];

profil = ProfilY(bb1, M);
tendance = filter(b, a, profil); % tendance décalée à cause du retard de groupe

[phi, wPhi] = phasez(b, a, M); 
retardPhi = abs(floor((phi(20) - phi(2)) / (wPhi(20) - wPhi(2)))); % dérivée de la phase
tendanceRecalee = tendance(retardPhi:size(profil,2));

%plot profil seul
% subplot(3,1,1), plot(profil), xlim([0, 2001]);
% xlim([0;M]);
% xticks([0,500,1000,1500,2000])
% xticklabels({'0', '0.5', '1', '1.5', '2'});
% ylabel("A");
% xlabel("Temps (s)");

%plot profil + fit décalé
subplot(2,1,1), xlim([0, 2001]);
hold on
plot(profil);
plot(tendance, 'r');
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});
ylabel("A");
xlabel("Temps (s)");
hold off

%plot profil + fit recalé
subplot(2,1,2), xlim([0, 2001]);
hold on
plot(profil);
plot(tendanceRecalee, 'r');
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});
ylabel("A");
xlabel("Temps (s)");
hold off

%% représentations du filtre 

clc
clear all
close all
warning('off', 'all');

M = 2001;
N = 27;
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
xlabel('Pulsation w (\times\pi rad/sample)')
ylabel('Gain (dB)')

[phi, wPhi] = phasez(b, a, M);
figure, plot(wPhi, phi);
ylabel("Phase (rad)");
xlabel("Pulsation w (\times\pi rad/sample)");
hold off

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

hold on;
plot(logN1, logF1, 'or');
plot(logN1, fit1);
ylabel("log(F)");
xlabel("log(N)");
title("Représentation des log de F par DMA");
text(logN1(9), logF1(9), sprintf("   alpha = %.3f", alpha1));
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
ylabel("Alpha");
title("Alphas moyens par DMA");

%% statistiques DFA/DMA

clc
clear all
warning('off', 'all');
close all

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
ylabel("Alpha");
title("Alphas moyens");

%%  TF + spectrogramme signal réel

clc
clear
close all

M = 2001;
data = load('dataEEG2020.mat');
y = cell2mat(data.dataEEG2020e7(2,1))';


% plot signal 
figure,
subplot(3,1,1);
plot(y);
xlabel('Temps en secondes');
ylabel('Intensité signal');
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});

% TF
tf = fft(y);

subplot(3,1,2);
plot(fftshift(abs(tf)));
xlabel('Fréquences normalisées');
ylabel('|TF(y)|');
xlim([0;M]);
xticks([0,1000,2000])
xticklabels({'-0.5', '0', '0.5'});

% spectrogram
subplot(3,1,3);
spectrogram(y - mean(y), sqrt(M), 'yaxis');
xlabel('Temps en secondes');
ylabel('Frequences');
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});

%% comparaison spectrogrammes entre deux états

clc
clear
close all

M = 2001;
data = load('dataEEG2020.mat');
y1 = cell2mat(data.dataEEG2020e7(1,1))';
y2 = cell2mat(data.dataEEG2020e7(2,1))';

figure,
subplot(2,1,1);
spectrogram(y1 - mean(y1), sqrt(M), 'yaxis');
xlabel('Temps en secondes');
ylabel('Fréquences');
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});
title('Etat 1');

subplot(2,1,2);
spectrogram(y2 - mean(y2), sqrt(M), 'yaxis');
xlabel('Temps en secondes');
ylabel('Frequences');
xlim([0;M]);
xticks([0,500,1000,1500,2000])
xticklabels({'0', '0.5', '1', '1.5', '2'});
title('Etat 2');

%% analyse signaux réels

clc
clear
close all

data = load('dataEEG2020.mat');

%création du dataframe
alpha = zeros(28, 1);
electrode = repmat(" ", 28,1);
method = repmat("   ", 28,1);
etat = repmat(' ', 28,1);
df = table(electrode, method, etat, alpha);
% sujet = repmat(' ', 28,1);

% df = table(sujet, electrode, etat, alpha);
iTableau = 1;


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
        a = Regularite(y, 'DFA');
        alphaDFAe1(1, index) = a;
        index = index +1;
        
        if(p == 1)
            tempEtat = 'A';
        else
            tempEtat = 'B';
        end
        
        df.electrode(iTableau) = "e1";
        df.method(iTableau) = {'DFA'};
        df.etat(iTableau) = tempEtat;
        df.alpha(iTableau) = a;
        iTableau = iTableau + 1;
    end
    
end

% DFA electrode 2
alphaDFAe2 = zeros(1, nbPhases*nbSujets);
index = 1;
for p = 1:nbPhases
    for s = 1:nbSujets
        
        % le signal est à transposer car on l'obtient en vecteur colonne
        y = cell2mat(data.dataEEG2020e8(p,s))'; 
        a = Regularite(y, 'DFA');
        alphaDFAe2(1, index) = a;
        index = index +1;
        
        if(p == 1)
            tempEtat = 'A';
        else
            tempEtat = 'B';
        end
        
        df.electrode(iTableau) = "e2";
        df.method(iTableau) = 'DFA';
        df.etat(iTableau) = tempEtat;
        df.alpha(iTableau) = a;
        iTableau = iTableau + 1;
    end
    
end



% analyse avec DMA -------------------------------

nbElectrodes = 2;
nbPhases = size(data.dataEEG2020e7, 1);
nbSujets = size(data.dataEEG2020e7, 2);

% DMA electrode 1
alphaDMAe1 = zeros(1, nbPhases*nbSujets);
index = 1;
for p = 1:nbPhases
    for s = 1:nbSujets
        
        y = cell2mat(data.dataEEG2020e7(p,s))';
        a = Regularite(y, 'DMA');
        alphaDMAe1(1, index) = a;
        index = index +1;
        
        if(p == 1)
            tempEtat = 'A';
        else
            tempEtat = 'B';
        end
        
        df.electrode(iTableau) = "e1";
        df.method(iTableau) = 'DMA';
        df.etat(iTableau) = tempEtat;
        df.alpha(iTableau) = a;
        iTableau = iTableau + 1;
    end
    
end

% DMA electrode 2
alphaDMAe2 = zeros(1, nbPhases*nbSujets);
index = 1;
for p = 1:nbPhases
    for s = 1:nbSujets
        
        % le signal est à transposer car on l'obtient en vecteur colonne
        y = cell2mat(data.dataEEG2020e8(p,s))'; 
        a = Regularite(y, 'DMA');
        alphaDMAe2(1, index) = a;
        index = index +1;
        
        if(p == 1)
            tempEtat = 'A';
        else
            tempEtat = 'B';
        end
        
        df.electrode(iTableau) = "e2";
        df.method(iTableau) = 'DMA';
        df.etat(iTableau) = tempEtat;
        df.alpha(iTableau) = a;
        iTableau = iTableau + 1;
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



display(df);
% writetable(df, "data_eeg.csv");

%% scatterplots

close all

%e1 DFA
a = df(df.electrode=="e1" & df.etat=='A' & df.method=="DFA", :).alpha;
b = df(df.electrode=="e1" & df.etat=='B' & df.method=="DFA", :).alpha;

figure,
subplot(2,1,1);
hold on
    plot(a, "r");
    plot(b, "g");
    yline(0.85, 'p');
    xlim([0,8]);
    ylabel("Alphas");
    title("Regularité DFA e1");
    legend("e1 etat A", "e1 etat B");
hold off

%e2 DFA
a = df(df.electrode=="e2" & df.etat=='A' & df.method=="DFA", :).alpha;
b = df(df.electrode=="e2" & df.etat=='B' & df.method=="DFA", :).alpha;

subplot(2,1,2);
hold on
    plot(a, "r");
    plot(b, "g");
    xlim([0,8]);
    ylabel("Alphas");
    title("Regularité DFA e2");
    legend("e2 etat A", "e2 etat B");
hold off

%e1 DMA
a = df(df.electrode=="e1" & df.etat=='A' & df.method=="DMA", :).alpha;
b = df(df.electrode=="e1" & df.etat=='B' & df.method=="DMA", :).alpha;

figure,
subplot(2,1,1);
hold on
    plot(a, "r");
    plot(b, "g");
    xlim([0,8]);
    ylabel("Alphas");
    title("Regularité DMA e1");
    legend("e1 etat A", "e1 etat B");
hold off

%e2 DMA
a = df(df.electrode=="e2" & df.etat=='A' & df.method=="DMA", :).alpha;
b = df(df.electrode=="e2" & df.etat=='B' & df.method=="DMA", :).alpha;

subplot(2,1,2);
hold on
    plot(a, "r");
    plot(b, "g");
    xlim([0,8]);
    ylabel("Alphas");
    title("Regularité DMA e2");
    legend("e2 etat A", "e2 etat B");
hold off

