function out = conv_fft(x, k, shape)
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
