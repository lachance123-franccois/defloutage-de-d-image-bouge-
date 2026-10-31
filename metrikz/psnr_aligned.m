function [p, shift] = psnr_aligned(ref, est, m)
%PSNR_ALIGNED PSNR après recalage entier (décalage de +-m pixels au plus).
%   En déconvolution aveugle, le noyau n'est connu qu'à une translation
%   près : l'image restaurée peut être décalée de quelques pixels. On
%   évalue donc au meilleur décalage, en ignorant une marge de m pixels
%   (pratique standard, Levin et al. 2009).
    if nargin < 3, m = 8; end
    p = -inf; shift = [0 0];
    a = ref(m+1:end-m, m+1:end-m);
    for dy = -m:m
        for dx = -m:m
            b = circshift(est, [dy dx]);
            v = psnr_db(a, min(max(b(m+1:end-m, m+1:end-m), 0), 1));
            if v > p
                p = v; shift = [dy dx];
            end
        end
    end
end
