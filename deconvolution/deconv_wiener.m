function x = deconv_wiener(y, k, nsr)
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
