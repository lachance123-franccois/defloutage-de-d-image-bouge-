function k = shake_kernel(n)
    if nargin >= 2
        rng(100);
    end
    T = 3000;
    t = linspace(0, 1, T);
    px = zeros(1, T);
    py = zeros(1, T);
    for j = 1:3
        f = 0.2 + rand;                                  
        px = px + randn * sin(2 * pi * f * t + 2 * pi * rand);
        py = py + randn * sin(2 * pi * f * t + 2 * pi * rand);
    end
    ext = max(max(px) - min(px), max(py) - min(py));
    s = 0.75 * (n - 3) / ext;
    px = (px - (min(px) + max(px)) / 2) * s;
    py = (py - (min(py) + max(py)) / 2) * s;
    c = (n + 1) / 2;
    k = splat_points(zeros(n, n), c + px, c + py);
    k = k / sum(k(:));
end
