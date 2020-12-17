function [F] = CalculFDMA( signal, M, N )
% Calcul le le F d'un signal grâce à une DMA

profil = ProfilY(signal, M);

b = zeros(1, N+1);
b(1) = 1;
b(N+1) = -1;
a = [N, -N];
tendance = filter(b, a, profil); % tendance décalée à cause du retard de groupe

[phi, wPhi] = phasez(b, a, M); 
retardPhi = abs(round((phi(20) - phi(2)) / (wPhi(20) - wPhi(2)))); % dérivée de la phase
tendanceRecalee = tendance(retardPhi:size(profil,2));
residus = profil(1:size(tendanceRecalee, 2)) - tendanceRecalee;

F = sqrt((1 / M)*sum(residus .^2));
end

