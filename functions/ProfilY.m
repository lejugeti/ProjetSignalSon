function [profil] = ProfilY(signal, M)

mu = mean(signal);
profil = zeros(1, M);

for m = 1:M
    tempY = signal(1:m);
    profil(1, m) = sum(tempY - mu);
    
end

end

