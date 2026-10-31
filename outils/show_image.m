function show_image(I, ttl, lims)
%SHOW_IMAGE Affiche une image en niveaux de gris ou RGB, sans toolbox.
%   show_image(I)               I en double dans [0, 1]
%   show_image(I, titre)
%   show_image(I, titre, [vmin vmax])   plage d'affichage
    if nargin < 3 || isempty(lims)
        lims = [0 1];
    end
    I = double(I);
    if size(I, 3) == 3
        image(min(max((I - lims(1)) / (lims(2) - lims(1)), 0), 1));
    else
        imagesc(I, lims);
        colormap(gca, gray(256));
    end
    axis image off;
    if nargin >= 2 && ~isempty(ttl)
        title(ttl, 'Interpreter', 'none', 'FontSize', 9);
    end
end
