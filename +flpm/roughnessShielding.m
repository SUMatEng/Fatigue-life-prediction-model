function [YR, sqrtAreaR] = roughnessShielding(Rv, RSm)
%FLPM.ROUGHNESSSHIELDING  Geometry factor of a rough surface (Murakami/Barricelli).
%
%   [YR, sqrtAreaR] = flpm.roughnessShielding(Rv, RSm)
%
%   Rv   mean valley depth a, RSm  mean pitch 2b of the roughness profile
%   (same length unit). The equivalent defect size of periodic surface
%   notches (Murakami) is
%     sqrt(areaR)/2b = 2.97 x - 3.51 x^2 - 9.74 x^3,  x = a/2b <= 0.195
%     sqrt(areaR)/2b = 0.38,                         0.195 < x <= 3
%   and the average geometry factor of the surface when the crack size is
%   taken as Rv is (Mater. Des. 2026, Eq. 9; Barricelli et al. 2023)
%     Y_R = 0.65 sqrt(sqrt(areaR) / a).
%
%   NOTE: Mater. Des. Eq. 7 prints the first coefficient as 0.97; the
%   published code (roughness_equivalent.m) and Murakami's original use 2.97,
%   which is adopted here (see docs/REVIEW.md). Vectorised.
x = Rv ./ RSm;
f = 0.38 * ones(size(x));
lowSlope = x > 0 & x <= 0.195;
f(lowSlope) = 2.97*x(lowSlope) - 3.51*x(lowSlope).^2 - 9.74*x(lowSlope).^3;
f(~(x > 0 & x <= 3)) = NaN;
sqrtAreaR = RSm .* f;
YR = 0.65 * sqrt(sqrtAreaR ./ Rv);
end
