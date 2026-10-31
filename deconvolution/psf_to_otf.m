function K = psf_to_otf(k, sz)
%PSF_TO_OTF Fonction de transfert d'un noyau centré (équivalent de psf2otf).
%   Le noyau est placé dans une image de taille sz puis décalé
%   circulairement pour que son centre soit en (1, 1).
    big = zeros(sz);
    [kh, kw] = size(k);
    big(1:kh, 1:kw) = k;
    big = circshift(big, [-floor(kh / 2), -floor(kw / 2)]);
    K = fft2(big);
end
