function out = resize_bilinear(x, outSize)
    [H, W] = size(x);
    r = ((1:outSize(1)) - 0.5) * H / outSize(1) + 0.5;
    c = ((1:outSize(2)) - 0.5) * W / outSize(2) + 0.5;
    r = min(max(r, 1), H);
    c = min(max(c, 1), W);
    [C, Rr] = meshgrid(c, r);
    out = interp2(x, C, Rr, 'linear');
end
