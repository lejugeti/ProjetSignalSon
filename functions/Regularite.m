function [alpha, beta, logN, logF] = Regularite( signal, method )
% calcule le coefficient de régularité d'un signal

M = size(signal, 2);
Ns = [13, 17, 21, 27, 35, 47, 59, 77, 99, 135, 189];
Fs = zeros(size(Ns));

if(method == 'DFA')
    for i = 1:size(Ns,2)
       N = Ns(i);
       Fs(i) = CalculF(signal, M, N);
    end

elseif(method == 'DMA')
      for i = 1:size(Ns,2)
          N = Ns(i);
          Fs(i) = CalculFDMA(signal, M, N);
      end
end

% calcul du fit pour obtenir le coefficient de régularité alpha
logF = log(Fs);
logN = log(Ns);
coefsDroite = polyfit(logN, logF, 1);

alpha = coefsDroite(1);
beta = coefsDroite(2);
end

