function J = edge_taper(I, k)
%EDGE_TAPER Prépare une image pour une déconvolution par FFT (cf. edgetaper).
%   La FFT suppose l'image périodique : la discontinuité entre bords opposés
%   crée de fortes oscillations (ringing). On remplace progressivement,
%   près des bords, l'image par sa version floutée circulairement, qui est
%   cohérente avec le modèle circulaire.
    [H, W] = size(I);
    kh = size(k, 1);
    w = ramp(H, kh) * ramp(W, kh)';
    Ib = real(ifft2(fft2(I) .* psf_to_otf(k, [H W])));
    J = w .* I + (1 - w) .* Ib;
end

function w = ramp(n, kh)
    w = ones(n, 1);
    m = min(kh, floor(n / 4));
    t = ((0:m-1)' + 0.5) / m;
    w(1:m) = 0.5 - 0.5 * cos(pi * t);
    w(end-m+1:end) = flipud(w(1:m));
end
