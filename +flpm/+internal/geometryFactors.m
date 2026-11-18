function [Y, cls] = geometryFactors(S, D, opts, rScale)
%GEOMETRYFACTORS  Stress-intensity geometry factor Y of every defect.
%   [Y, cls] = flpm.internal.geometryFactors(S, D, opts, rScale)
%   D in m. cls: 0 = internal, 1 = near-surface, 2 = surface-breaking.
%   Classification (Murakami proximity rule, discretizeflaws.m):
%     surface       Rpos + r >= D/2
%     near-surface  r / (D/2 - Rpos - r) >= opts.NearSurfaceRatio
%     internal      otherwise
%   with r = rScale*a. Embedded defects use the Alegre et al. (2022) solution
%   for an embedded elliptical crack (the larger vertex value); internal ones
%   are further scaled by opts.InternalYFactor (0.5/0.65, Murakami). Surface-
%   breaking defects use the Forman-Shivakumar surface-crack solution,
%   bounded below by the defect's yFloor (roughness shielding factor).
%   The Alegre 'h' argument is D/2 - Rpos + a, as in the published code.
Rb = D/2;
Rpos = hypot(S.x, S.y);
r = rScale * S.a;
Rr = Rpos + r;
surf = Rr >= Rb;
lig = Rb - Rr;
near = ~surf & (r ./ lig >= opts.NearSurfaceRatio);
Y = zeros(size(S.a));
emb = ~surf;
if any(emb)
    Y(emb) = flpm.internal.alegreEmbeddedY(S.a(emb), S.asp(emb), ...
        Rb - Rpos(emb) + S.a(emb), Rb);
end
internal = emb & ~near;
Y(internal) = opts.InternalYFactor * Y(internal);
if any(surf)
    Y(surf) = max(flpm.internal.formanShaftY(S.a(surf), D), S.yfloor(surf));
end
cls = zeros(size(S.a));
cls(near) = 1;
cls(surf) = 2;
end
