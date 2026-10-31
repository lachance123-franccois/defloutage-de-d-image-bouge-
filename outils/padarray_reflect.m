function P = padarray_reflect(f, r)
    [H, W, ~] = size(f);
    ri = [r+1:-1:2, 1:H, H-1:-1:H-r];
    ci = [r+1:-1:2, 1:W, W-1:-1:W-r];
    P = f(ri, ci, :);
end
