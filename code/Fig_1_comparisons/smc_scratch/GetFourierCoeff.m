function [muk, HK] = GetFourierCoeff(opt, X, Y, mu)
%% returns the fourier coeffients of the prior information distribution
Nkx = opt.erg.Nkx;
Nky = opt.erg.Nky;
Lx = opt.L(1);
Ly = opt.L(2);

% Initializing Fourier coefficents of prior
muk = zeros(Nkx,Nky);
HK = sqrt(Lx*Ly) .* [1 sqrt(1/2)*ones(1,Nky-1)]' * [1 sqrt(1/2)*ones(1,Nkx-1)];
mu = reshape(mu,size(X));

for kx = 0:Nkx-1
    for ky = 0:Nky-1
        muk(kx+1, ky+1) = sum(sum( mu .* cos(kx * pi * X/Lx) .* cos(ky * pi *  Y/Ly) )) / HK(ky+1,kx+1);
    end
end
end