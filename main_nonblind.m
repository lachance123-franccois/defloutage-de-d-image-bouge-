
clear; close all;
outDir = 'results';
if ~exist(outDir, 'dir'), mkdir(outDir); end
K = load(fullfile('data', 'kernels.mat'));
knames = {'K1', 'K2', 'K3', 'K4'};
kdesc = {'linéaire 15 px, 30°', 'linéaire 25 px, 120°', 'tremblement 1', 'tremblement 2'};
names = {'camera', 'astronaut', 'coffee', 'brick'};
methods = {'Flou', 'Wiener', 'Richardson-Lucy', 'TV'};
sigma_n = 0.01;
P = zeros(numel(knames), numel(names), numel(methods)); S = P;

for a = 1:numel(knames)
    k = K.(knames{a});
    for b = 1:numel(names)
        x = read_image(fullfile('data', [names{b} '.png']));
        rng(1);
        y = blur_image(x, k) + sigma_n * randn(size(x));
        outs = {y, deconv_wiener(y, k, 0.02), deconv_rl(y, k, 40), deconv_tv(y, k, 1e-3, 300)};
        for m = 1:numel(methods)
            u = min(max(outs{m}, 0), 1);
            P(a, b, m) = psnr_db(x, u);
            S(a, b, m) = ssim_index(x, u);
        end
        fprintf('%s %-10s', knames{a}, names{b});
        fprintf('  %s %.2f', methods{1}, P(a, b, 1));
        for m = 2:numel(methods)
            fprintf(' | %s %.2f dB / %.3f', methods{m}, P(a, b, m), S(a, b, m));
        end
        fprintf('\n');
        if b == 1                                   % figure sur cameraman
            fig = figure('Color', 'w', 'Position', [40 40 1400 330]);
            subplot(1, 6, 1); show_image(k / max(k(:)), sprintf('Noyau %s', knames{a}));
            subplot(1, 6, 2); show_image(x, 'Original');
            for m = 1:numel(methods)
                subplot(1, 6, 2 + m);
                show_image(outs{m}, sprintf('%s  %.2f dB', methods{m}, P(a, b, m)));
            end
            sgtitle(sprintf('Noyau %s (%s), bruit \\sigma = %.2f', knames{a}, kdesc{a}, sigma_n));
            save_figure(fig, fullfile(outDir, sprintf('nonblind_%s.png', knames{a})));
            close(fig);
        end
    end
end

fprintf('\nPSNR moyen sur les 4 images (dB)\n%-6s', '');
fprintf('%18s', methods{:}); fprintf('\n');
for a = 1:numel(knames)
    fprintf('%-6s', knames{a}); fprintf('%18.2f', squeeze(mean(P(a, :, :), 2))); fprintf('\n');
end

% ---------------- semi-convergence de Richardson-Lucy ----------------
x = read_image(fullfile('data', 'camera.png'));
k = K.K4;
rng(1);
y = blur_image(x, k) + sigma_n * randn(size(x));
its = [5 10 20 30 40 60 80 120 160];
pr = zeros(size(its));
for t = 1:numel(its)
    pr(t) = psnr_db(x, min(max(deconv_rl(y, k, its(t)), 0), 1));
end
lams = [2e-4 5e-4 1e-3 2e-3 5e-3 1e-2];
pt = zeros(size(lams));
for t = 1:numel(lams)
    pt(t) = psnr_db(x, deconv_tv(y, k, lams(t), 300));
end
fig = figure('Color', 'w', 'Position', [100 100 1000 380]);
subplot(1, 2, 1); plot(its, pr, '-o', 'LineWidth', 1.5); grid on;
xlabel('itérations'); ylabel('PSNR (dB)'); title('Richardson-Lucy : semi-convergence');
subplot(1, 2, 2); semilogx(lams, pt, '-o', 'LineWidth', 1.5); grid on;
xlabel('\lambda'); ylabel('PSNR (dB)'); title('TV : choix de \lambda');
save_figure(fig, fullfile(outDir, 'parameters.png'));
