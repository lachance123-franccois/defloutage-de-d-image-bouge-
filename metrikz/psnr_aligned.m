function [p, shift] = psnr_aligned(ref, est, m)
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
