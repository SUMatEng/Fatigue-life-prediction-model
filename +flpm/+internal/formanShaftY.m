function Y = formanShaftY(a, D)
%FORMANSHAFTY  Geometry factor of a surface crack in a round bar.
%   Y = flpm.internal.formanShaftY(a, D) with crack depth a and bar diameter
%   D in the same units. Forman & Shivakumar (ASTM STP 905, 1986) form as
%   implemented in the published code (surface_crack_shaft_geo.m): a is held
%   at D/2 for larger cracks. Vectorised; NaN in, NaN out.
a = min(a, D/2);
B = (pi/2) * (a ./ D);
G = 0.92 * (2/pi) .* sec(B) .* sqrt(tan(B) ./ B);
Z = 1 - sin(B);
Y = G .* (0.752 + 1.286 .* B + 0.37 .* Z.^3);
small = B < 1e-8;                     % limit a -> 0 (tan(B)/B -> 1)
Y(small) = 0.92*(2/pi)*(0.752 + 0.37);
end
