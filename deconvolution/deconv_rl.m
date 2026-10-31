function x = deconv_rl(y, k, n_iter)
    if nargin < 3
        n_iter = 40;
    end
    r = (size(k, 1) - 1) / 2;
    kf = rot90(k, 2);
    y = max(y, 1e-4);                                   
    [H, W] = size(y);
    x = y([ones(1, r), 1:H, H * ones(1, r)], [ones(1, r), 1:W, W * ones(1, r)]);
    nrm = conv_fft(ones(H, W), kf, 'full');            
    seen = nrm > 0.1;                                    
    nrm = max(nrm, 0.1);
    for it = 1:n_iter
        ratio = y ./ max(conv_fft(x, k, 'valid'), 1e-4);
        upd = min(x .* max(conv_fft(ratio, kf, 'full'), 0) ./ nrm, 1);
        x(seen) = upd(seen);
    end
    x = x(r+1:end-r, r+1:end-r);
end
