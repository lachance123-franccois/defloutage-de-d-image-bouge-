function x = deconv_wiener(y, k, nsr)
%DECONV_WIENER Déconvolution de Wiener (rapport bruit/signal constant).
%   x = deconv_wiener(y, k)         nsr = 0.02
%   x = deconv_wiener(y, k, nsr)
%
%   X(w) = conj(K(w)) Y(w) / (|K(w)|^2 + nsr)
%   nsr joue le rôle de régularisation : là où |K| est petit (zéros du
%   spectre du flou), on n'amplifie pas le bruit. Bords : extension miroir
%   puis edge_taper, pour limiter les oscillations dues à la FFT.
    if nargin < 3
        nsr = 0.02;
    end
    r = size(k, 1);
    yp = edge_taper(padarray_reflect(y, r), k);
    K = psf_to_otf(k, size(yp));
    X = conj(K) .* fft2(yp) ./ (abs(K).^2 + nsr);
    x = real(ifft2(X));
    x = x(r+1:end-r, r+1:end-r);
end
