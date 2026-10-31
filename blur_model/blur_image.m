function y = blur_image(x, k)
    r = (size(k, 1) - 1) / 2;
    [H, W, C] = size(x);
    ri = [ones(1, r), 1:H, H * ones(1, r)];
    ci = [ones(1, r), 1:W, W * ones(1, r)];
    y = zeros(H, W, C);
    for c = 1:C
        y(:, :, c) = conv2(x(ri, ci, c), k, 'valid');
    end
end
