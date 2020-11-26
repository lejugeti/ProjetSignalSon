function [F] = CalculProfilGlobal( signal, M, N )
% calcul du profil global F

L = floor(M / N);
profilY = ProfilY(signal, L*N);
coefsX = CalculCoefsProfil(profilY, L, N);
newY = zeros(L*N, 1); % vecteur des résidus

for l = 1:L
    for n = 1:N
        m = (l-1) * N + n;
        tempX = polyval(coefsX(l,:), m);
        newY(m, 1) = profilY(m) - tempX;
    end
end

F = sqrt( (1/(L*N))  * sum(newY .^ 2));
end