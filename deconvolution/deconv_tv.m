function x = deconv_tv(y, k, lambda, n_iter)
    if nargin < 3 || isempty(lambda)
        lambda = 1e-3;
    end
    if nargin < 4
        n_iter = 300;
    end
    r = (size(k, 1) - 1) / 2;
    kf = rot90(k, 2);
    [H, W] = size(y);
    x = y([ones(1, r), 1:H, H * ones(1, r)], [ones(1, r), 1:W, W * ones(1, r)]);
    x = min(max(x, 0), 1);
    xb = x;
    p = zeros(H, W);
    qx = zeros(size(x));
    qy = zeros(size(x));
    tau = 1 / 3;
    sigma = 1 / 3;
    for it = 1:n_iter
        % dual de l'attache aux données : prox de F*, F(z) = 1/2 ||z - y||^2
        p = (p + sigma * (conv_fft(xb, k, 'valid') - y)) / (1 + sigma);
        % dual de la TV : projection sur la boule ||q|| <= lambda
        [gx, gy] = grad_fwd(xb);
        qx = qx + sigma * gx;
        qy = qy + sigma * gy;
        nq = max(1, sqrt(qx.^2 + qy.^2) / lambda);
        qx = qx ./ nq;
        qy = qy ./ nq;
        % primal : pas de descente puis projection sur [0, 1]
        x_old = x;
        x = min(max(x - tau * (conv_fft(p, kf, 'full') - div_bwd(qx, qy)), 0), 1);
        xb = 2 * x - x_old;
    end
    x = x(r+1:end-r, r+1:end-r);
end

function [gx, gy] = grad_fwd(x)
%   gradient par différences avant (condition de Neumann)
    gx = [x(:, 2:end) - x(:, 1:end-1), zeros(size(x, 1), 1)];
    gy = [x(2:end, :) - x(1:end-1, :); zeros(1, size(x, 2))];
end

function d = div_bwd(px, py)
%   divergence = -adjoint de grad_fwd
    dx = [px(:, 1), px(:, 2:end-1) - px(:, 1:end-2), -px(:, end-1)];
    dy = [py(1, :); py(2:end-1, :) - py(1:end-2, :); -py(end-1, :)];
    d = dx + dy;
end
