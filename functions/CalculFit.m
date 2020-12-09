function [fitEntier] = CalculFit(signal, M, N, L)
% Calcule et renvoie le fit du profil du signal par des polynômes 
% de degré 3

profil = ProfilY(signal, M);
coefs = CalculCoefsProfil(profil, L, N);

fitEntier = zeros(1,L*N);
m = 1:N;
for l = 1:L
    
    tempCoefs = coefs(l, :); % coefs pour le l donné
    fit = polyval(tempCoefs, m); % calcul du fit pour tous les m
    fitEntier(m) = fit; % attribution du fit obtenu au fit entier
    m = m + N;
end
end

