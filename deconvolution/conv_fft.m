function out = conv_fft(x, k, shape)
%CONV_FFT Convolution 2D par FFT, même résultat que conv2(x, k, shape).
%   shape : 'full' ou 'valid'. Beaucoup plus rapide que conv2 pour les
%   grands noyaux (bougé de 25 pixels et plus).
    [H, W] = size(x);
    [kh, kw] = size(k);
    Hf = H + kh - 1;
    Wf = W + kw - 1;
    F = real(ifft2(fft2(x, Hf, Wf) .* fft2(k, Hf, Wf)));
    switch shape
        case 'full'
            out = F;
        case 'valid'
            out = F(kh:H, kw:W);
        otherwise
            error('conv_fft : forme ''%s'' non gérée', shape);
    end
end
