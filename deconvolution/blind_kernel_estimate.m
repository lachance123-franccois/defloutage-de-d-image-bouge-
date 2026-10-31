function k = blind_kernel_estimate(y, ks, lambda, gamma, n_inner)s
    if nargin < 3, lambda = 4e-3; end
    if nargin < 4, gamma = 100; end
    if nargin < 5, n_inner = 10; end
    n_scales = max(1, floor(log(ks / 5) / log(sqrt(2))) + 1);
    k = [];
    for s = n_scales-1:-1:0
        f = sqrt(0.5)^s;
        kss = ceil(ks * f);
        if mod(kss, 2) == 0, kss = kss + 1; end
        if s > 0
            ys = resize_bilinear(y, round(size(y) * f));
        else
            ys = y;
        end
        if isempty(k)
            k = zeros(kss);                    
            c = (kss + 1) / 2;
            k(c, c-1:c+1) = 1 / 3;
            k(c-1:c+1, c) = k(c-1:c+1, c) + 1 / 3;
            k = k / sum(k(:));
        else
            k = clean_kernel(resize_bilinear(k, [kss kss]));
        end
        pad = kss;
        yp = edge_taper(padarray_reflect(ys, pad), k);
        m = zeros(size(yp));
        m(pad+1:end-pad, pad+1:end-pad) = 1;   % on ignore la marge artificielle
        [gyx, gyy] = grad_circ(yp);
        lam = lambda;
        for it = 1:n_inner
            x = l0_restore(yp, k, lam);
            [gxx, gxy] = grad_circ(x);
            k = estimate_psf({gyx .* m, gyy .* m}, {gxx .* m, gxy .* m}, kss, gamma);
            k = center_kernel(clean_kernel(k));
            lam = max(lam / 1.1, 1e-4);
        end
    end
end

% -------------------------------------------------------------------------
function [gx, gy] = grad_circ(x)
    gx = circshift(x, [0 -1]) - x;
    gy = circshift(x, [-1 0]) - x;
end

function x = l0_restore(y, k, lambda)
%   min_x ||k*x - y||^2 + lambda ||grad x||_0  (Xu et al. 2011) :
%   décomposition quadratique avec beta croissant ; à beta fixé, le
%   gradient auxiliaire est seuillé (L0), puis x est obtenu par FFT.
    sz = size(y);
    K = psf_to_otf(k, sz);
    Dx = psf_to_otf([1 -1], sz);
    Dy = psf_to_otf([1; -1], sz);
    KtY = conj(K) .* fft2(y);
    KtK = abs(K).^2;
    DtD = abs(Dx).^2 + abs(Dy).^2;
    x = y;
    beta = 2 * lambda;
    while beta < 1e5
        [gx, gy] = grad_circ(x);
        t = (gx.^2 + gy.^2) < lambda / beta;
        gx(t) = 0;
        gy(t) = 0;
        dv = (circshift(gx, [0 1]) - gx) + (circshift(gy, [1 0]) - gy);   % adjoint du gradient
        x = real(ifft2((KtY + beta * fft2(dv)) ./ (KtK + beta * DtD)));
        beta = beta * 2;
    end
end

function k = estimate_psf(gy, gx, ks, gamma)
%   min_k  somme_d ||k * gx_d - gy_d||^2 + gamma ||k||^2   (Fourier)
    num = 0; den = 0;
    for d = 1:numel(gy)
        BX = fft2(gx{d});
        BY = fft2(gy{d});
        num = num + conj(BX) .* BY;
        den = den + abs(BX).^2;
    end
    kf = fftshift(real(ifft2(num ./ (den + gamma))));
    [H, W] = size(kf);
    cy = floor(H / 2) + 1;
    cx = floor(W / 2) + 1;
    r = (ks - 1) / 2;
    k = kf(cy-r:cy+r, cx-r:cx+r);
end

function k = clean_kernel(k)
%   positivité, seuil à 5 % du max, composantes connexes de masse < 10 %
%   supprimées (bruit), normalisation
    k = max(k, 0);
    if sum(k(:)) == 0, return; end
    k = k / sum(k(:));
    k(k < 0.05 * max(k(:))) = 0;
    lab = label8(k > 0);
    for c = 1:max(lab(:))
        if sum(k(lab == c)) < 0.1
            k(lab == c) = 0;
        end
    end
    k = k / sum(k(:));
end

function lab = label8(bw)
%   étiquetage des composantes 8-connexes (parcours en largeur)
    [H, W] = size(bw);
    lab = zeros(H, W);
    n = 0;
    for p = find(bw)'
        if lab(p) > 0, continue; end
        n = n + 1;
        queue = p;
        lab(p) = n;
        while ~isempty(queue)
            q = queue(end);
            queue(end) = [];
            [i, j] = ind2sub([H W], q);
            for di = -1:1
                for dj = -1:1
                    a = i + di; b = j + dj;
                    if a >= 1 && a <= H && b >= 1 && b <= W && bw(a, b) && lab(a, b) == 0
                        lab(a, b) = n;
                        queue(end+1) = a + (b - 1) * H; %#ok<AGROW>
                    end
                end
            end
        end
    end
end

function k = center_kernel(k)
%   recentre le noyau sur son barycentre (décalage entier)
    n = size(k, 1);
    c = (n + 1) / 2;
    [xx, yy] = meshgrid(1:n, 1:n);
    cy = round(sum(yy(:) .* k(:)) / sum(k(:)));
    cx = round(sum(xx(:) .* k(:)) / sum(k(:)));
    k = circshift(k, [c - cy, c - cx]);
end
