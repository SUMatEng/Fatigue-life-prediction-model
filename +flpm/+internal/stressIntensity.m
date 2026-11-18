function [Y, A, cls] = stressIntensity(S)
%STRESSINTENSITY  Murakami sqrt(area) geometry of every defect (half-space).
%   [Y, A, cls] = flpm.internal.stressIntensity(S)
%   K = Y * sigma * sqrt(pi * A) with A the effective sqrt(area) (m).
%
%   Each defect is a disc of area S.a^2 (radius r = S.a/sqrt(pi)) in the
%   plane normal to the load, centred at depth S.z below the free surface.
%     surface-breaking  z <= r      A = sqrt(area of the disc below the
%                                       surface), Y = 0.65
%     near-surface      z - r <= r/0.8  A = sqrt(pi r (z + r)/2), the
%                                       half-ellipse enclosing the defect
%                                       and its ligament, Y = 0.65
%     internal          otherwise   A = S.a, Y = 0.5
%   No specimen dimension enters. cls: 0 internal, 1 near, 2 surface.
k = flpm.internal.constants();
r = S.a / sqrt(pi);
z = S.z;
surf = z <= r;
near = ~surf & (z - r) <= r / k.proximity;

A = S.a;
Y = k.Yinternal * ones(size(A));
if any(surf)
    rs = r(surf); zs = max(z(surf), -rs);
    inside = rs.^2 .* acos(-zs ./ rs) + zs .* sqrt(max(rs.^2 - zs.^2, 0));
    A(surf) = sqrt(inside);
    Y(surf) = k.Ysurface;
end
if any(near)
    A(near) = sqrt(pi * r(near) .* (z(near) + r(near)) / 2);
    Y(near) = k.Ysurface;
end
cls = zeros(size(A));
cls(near) = 1;
cls(surf) = 2;
end
