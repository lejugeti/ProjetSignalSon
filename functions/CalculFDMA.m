function [F] = CalculFDMA( signal, M, N )
% Calcul le le F d'un signal grâce à une DMA

profil = ProfilY(signal, M);
filtre = ones(1, N) / N;
tendance = conv(profil, filtre, 'same');
residus = profil - tendance;

F = (1 / M) * sqrt(sum(residus .^2));
end

