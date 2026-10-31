function x = deconv_rl(y, k, n_iter)
%DECONV_RL Richardson-Lucy avec gestion correcte des bords.
%   x = deconv_rl(y, k)            40 itérations
%   x = deconv_rl(y, k, n_iter)
%
%   Maximum de vraisemblance pour un bruit de Poisson, itérations
%   multiplicatives :   x <- x . A'( y ./ (A x) ) ./ A'(1)
%   A est la convolution 'valid' : x est plus grand que y d'une demi-taille
%   de noyau de chaque côté (la scène déborde du cadre), ce qui évite les
%   artefacts de bord (Bertero & Boccacci, 2005). A' = corrélation 'full'.
%   Le nombre d'itérations sert de régularisation : au-delà d'un certain
%   point, RL amplifie le bruit (semi-convergence).
    if nargin < 3
        n_iter = 40;
    end
    r = (size(k, 1) - 1) / 2;
    kf = rot90(k, 2);
    y = max(y, 1e-4);                                   % RL suppose des données positives
    [H, W] = size(y);
    x = y([ones(1, r), 1:H, H * ones(1, r)], [ones(1, r), 1:W, W * ones(1, r)]);
    nrm = conv_fft(ones(H, W), kf, 'full');             % 1 au centre, < 1 sur la marge
    seen = nrm > 0.1;                                    % pixels de marge assez observés
    nrm = max(nrm, 0.1);
    for it = 1:n_iter
        ratio = y ./ max(conv_fft(x, k, 'valid'), 1e-4);
        upd = min(x .* max(conv_fft(ratio, kf, 'full'), 0) ./ nrm, 1);
        x(seen) = upd(seen);
    end
    x = x(r+1:end-r, r+1:end-r);
end
