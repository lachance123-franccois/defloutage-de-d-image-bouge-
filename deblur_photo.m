clear; close all; rng(0);
photoFile = fullfile('data', 'camera.png');
mode      = 'aveugle';     
ks        = 31;            
maxSide   = 600;          
lambda    = 2e-3;       


I = read_image(photoFile);
s = max(size(I, 1), size(I, 2)) / maxSide;
if s > 1
    I = I(round(1:s:end), round(1:s:end), :);  
end

Y = 0.299 * I(:, :, 1) + 0.587 * I(:, :, 2) + 0.114 * I(:, :, 3);

tic;
switch mode
    case 'lineaire'
        [L, th] = estimate_motion_cepstrum(Y);
        k = motion_kernel(L, th);
        fprintf('Bougé estimé : %.1f px à %.0f°\n', L, th);
    otherwise
        k = blind_kernel_estimate(Y, ks);
end
fprintf('Noyau estimé en %.0f s\n', toc);

out = zeros(size(I));
for c = 1:size(I, 3)
    out(:, :, c) = deconv_tv(I(:, :, c), k, lambda, 300);
end
fig = figure('Color', 'w', 'Position', [50 50 1300 450]);
subplot(1, 3, 1); show_image(min(max(I, 0), 1), 'Photo floue');
subplot(1, 3, 2); show_image(k / max(k(:)), 'Noyau estimé');
subplot(1, 3, 3); show_image(out, 'Restaurée (TV)');
save_figure(fig, fullfile('results', 'my_photo.png'));
