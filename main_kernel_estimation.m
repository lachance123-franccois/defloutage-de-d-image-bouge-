clear; close all;
outDir = 'results';
if ~exist(outDir, 'dir'), mkdir(outDir); end
K = load(fullfile('data', 'kernels.mat'));
names = {'camera', 'astronaut', 'coffee', 'brick'};
sigma_n = 0.01;

% ---------------- 1. cepstre ----------------
truth = [15 30; 25 120];
fprintf('Bougé linéaire, estimation par le cepstre :\n');
for a = 1:2
    k = K.(sprintf('K%d', a));
    for b = 1:numel(names)
        x = read_image(fullfile('data', [names{b} '.png']));
        rng(1);
        y = blur_image(x, k) + sigma_n * randn(size(x));
        [L, th] = estimate_motion_cepstrum(y);
        ke = motion_kernel(L, th);
        u = deconv_tv(y, ke, 1e-3, 300);
        fprintf('  K%d (L = %d, theta = %d°) %-10s : L = %5.1f, theta = %5.1f°  -> TV %.2f dB\n', ...
            a, truth(a, 1), truth(a, 2), names{b}, L, th, psnr_aligned(x, u));
        if a == 2 && b == 1
            Yl = log(1 + abs(fftshift(fft2(y - mean(y(:))))));
            win = (0.5 - 0.5 * cos(2 * pi * (0:255)' / 255)) * (0.5 - 0.5 * cos(2 * pi * (0:255) / 255));
            Cp = fftshift(real(ifft2(log(1e-3 + abs(fft2((y - mean(y(:))) .* win))))));
            fig = figure('Color', 'w', 'Position', [100 100 1200 380]);
            subplot(1, 3, 1); show_image(y, 'Image floue (K2)');
            subplot(1, 3, 2); imagesc(Yl); axis image off; colormap(gca, gray);
            title('log |spectre| : bandes perpendiculaires au mouvement');
            subplot(1, 3, 3); imagesc(-Cp(129-40:129+40, 129-40:129+40), [0 max(-Cp(:)) / 4]);
            axis image off; colormap(gca, hot);
            title(sprintf('cepstre (inversé) : pics à L = %.1f, %.0f°', L, th));
            save_figure(fig, fullfile(outDir, 'cepstrum.png'));
            close(fig);
        end
    end
end

% ---------------- 2. estimation aveugle ----------------
fprintf('\nEstimation aveugle (prior L0, grossier -> fin) :\n');
fig = figure('Color', 'w', 'Position', [30 30 1300 1000]);
for a = 1:4
    k = K.(sprintf('K%d', a));
    sims = zeros(1, numel(names)); pb = sims; pt = sims;
    for b = 1:numel(names)
        x = read_image(fullfile('data', [names{b} '.png']));
        rng(1);
        y = blur_image(x, k) + sigma_n * randn(size(x));
        tic;
        ke = blind_kernel_estimate(y, 27);
        t = toc;
        sims(b) = kernel_similarity(ke, k);
        u = deconv_tv(y, ke, 2e-3, 300);
        pb(b) = psnr_aligned(x, u);
        pt(b) = psnr_db(x, deconv_tv(y, k, 1e-3, 300));
        fprintf('  K%d %-10s : similarité %.2f | TV noyau estimé %.2f dB | TV vrai noyau %.2f dB | flou %.2f dB | %.0f s\n', ...
            a, names{b}, sims(b), pb(b), pt(b), psnr_db(x, min(max(y, 0), 1)), t);
        if b == 1
            subplot(4, 4, 4 * (a - 1) + 1); show_image(k / max(k(:)), sprintf('K%d vrai', a));
            subplot(4, 4, 4 * (a - 1) + 2); show_image(ke / max(ke(:)), sprintf('estimé (sim. %.2f)', sims(b)));
            subplot(4, 4, 4 * (a - 1) + 3); show_image(y, 'flou');
            subplot(4, 4, 4 * (a - 1) + 4); show_image(u, sprintf('restauré %.2f dB', pb(b)));
        end
    end
end
sgtitle('Déconvolution aveugle : noyau estimé à partir de l''image floue seule');
save_figure(fig, fullfile(outDir, 'blind.png'));
