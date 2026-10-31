function J = edge_taper(I, k)
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
