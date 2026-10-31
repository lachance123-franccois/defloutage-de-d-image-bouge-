function s = kernel_similarity(k, kt)
    c = conv2(k, rot90(kt, 2), 'full');
    s = max(c(:)) / sqrt(sum(k(:).^2) * sum(kt(:).^2));
end
