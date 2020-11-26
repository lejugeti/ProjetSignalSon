function [coefsX] = CalculCoefsProfil( profilY, L, N )
% calcul des coefficients des polynômes qui fit le profil
% du signal

coefsX = zeros(L, 4);
for l = 1:L
    
    n = 1:N;
    m = (l - 1) * N + n; % vecteur
    tempY = profilY(m);
    coefs = polyfit(m, tempY, 3);
    coefsX(l,:) = coefs;
   
end

end

