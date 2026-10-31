function k = splat_points(k, xs, ys)
    x0 = floor(xs(:)); y0 = floor(ys(:));
    fx = xs(:) - x0;  fy = ys(:) - y0;
    [n1, n2] = size(k);
    ok = x0 >= 1 & x0 < n2 & y0 >= 1 & y0 < n1;
    x0 = x0(ok); y0 = y0(ok); fx = fx(ok); fy = fy(ok);
    k = k + accumarray([y0,     x0    ], (1 - fx) .* (1 - fy), [n1 n2]) ...
          + accumarray([y0,     x0 + 1], fx .* (1 - fy),       [n1 n2]) ...
          + accumarray([y0 + 1, x0    ], (1 - fx) .* fy,       [n1 n2]) ...
          + accumarray([y0 + 1, x0 + 1], fx .* fy,             [n1 n2]);
end
