function s = kernel_similarity(k, kt)
%KERNEL_SIMILARITY Corrélation normalisée maximale entre deux noyaux (1 = identiques
%   à une translation près, 0 = sans rapport).
    c = conv2(k, rot90(kt, 2), 'full');
    s = max(c(:)) / sqrt(sum(k(:).^2) * sum(kt(:).^2));
end
