function k = motion_kernel(L, theta)
    r = ceil(L / 2) + 1;
    n = 2 * r + 1;
    t = linspace(-L / 2, L / 2, ceil(10 * L) + 1);
    xs = r + 1 + t * cosd(theta);         
    ys = r + 1 - t * sind(theta);          
    k = splat_points(zeros(n, n), xs, ys);
    k = k / sum(k(:));
end
