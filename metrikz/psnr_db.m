function p = psnr_db(ref, est, peak)
    if nargin < 3
        peak = 1;
    end
    mse = mean((double(ref(:)) - double(est(:))).^2);
    p = 10 * log10(peak^2 / mse);
end
