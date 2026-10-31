function [L, theta] = estimate_motion_cepstrum(y, Lmax)
%ESTIMATE_MOTION_CEPSTRUM Estime longueur et angle d'un bougé linéaire.
%   [L, theta] = estimate_motion_cepstrum(y)
%   [L, theta] = estimate_motion_cepstrum(y, Lmax)     (défaut Lmax = 40)
%
%   Le spectre d'un flou linéaire de longueur L est un sinus cardinal : il
%   s'annule périodiquement (bandes sombres perpendiculaires au mouvement).
%   Le cepstre  C = ifft2( log|Y| )  transforme cette périodicité en pics
%   NÉGATIFS situés à +-L pixels de l'origine, dans la direction du
%   mouvement. On cherche le minimum de C dans une couronne 3 <= r <= Lmax.
%   (Rom 1975 ; Moghaddam & Jamzad 2007)
    if nargin < 2
        Lmax = 40;
    end
    [H, W] = size(y);
    win = hann_win(H) * hann_win(W)';                  % atténue les bords
    Y = abs(fft2((y - mean(y(:))) .* win));
    C = fftshift(real(ifft2(log(1e-3 + Y))));
    cy = floor(H / 2) + 1;                              % origine après fftshift
    cx = floor(W / 2) + 1;
    [xx, yy] = meshgrid(1:W, 1:H);
    d = hypot(yy - cy, xx - cx);
    C(d < 3 | d > Lmax) = inf;
    [~, idx] = min(C(:));
    [i, j] = ind2sub([H W], idx);
    dy = i - cy;
    dx = j - cx;
    L = hypot(dx, dy);
    theta = mod(atan2d(-dy, dx), 180);
end

function w = hann_win(n)
    w = 0.5 - 0.5 * cos(2 * pi * (0:n-1)' / (n - 1));
end
