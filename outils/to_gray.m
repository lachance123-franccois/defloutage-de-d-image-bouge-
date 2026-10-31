function g = to_gray(I).
    I = double(I);
    if size(I, 3) == 3
        g = 0.2989 * I(:, :, 1) + 0.5870 * I(:, :, 2) + 0.1140 * I(:, :, 3);
    else
        g = I;
    end
end
