function K = psf_to_otf(k, sz)
    big = zeros(sz);
    [kh, kw] = size(k);
    big(1:kh, 1:kw) = k;
    big = circshift(big, [-floor(kh / 2), -floor(kw / 2)]);
    K = fft2(big);
end
