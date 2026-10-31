function P = padarray_reflect(f, r)
%PADARRAY_REFLECT Extension miroir de r pixels, sans répéter le bord.
%   Exemple pour r = 2 : [c b | a b c d | c b].  Fonctionne aussi sur des
%   images multicanaux (r doit être < taille de l'image).
    [H, W, ~] = size(f);
    ri = [r+1:-1:2, 1:H, H-1:-1:H-r];
    ci = [r+1:-1:2, 1:W, W-1:-1:W-r];
    P = f(ri, ci, :);
end
