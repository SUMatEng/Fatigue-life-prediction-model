function s = pennyStress(rho, z)
%PENNYSTRESS  Normal stress of a penny-shaped crack under unit pressure.
%   s = pennyStress(rho, z) gives sigma_zz / p at in-plane distance rho
%   and normal distance z from the centre of a circular crack of unit
%   radius whose faces carry a uniform opening pressure p, in an infinite
%   isotropic body (Sneddon's solution, closed form):
%     sigma_zz/p = -(2/pi) [G - z dG/dz]
%     G     = asin(1/l2) - Re Q,           Q = (rho^2 + (z + i)^2)^(-1/2)
%     dG/dz = Im Q + Re[(z + i) Q^3]
%   with l2 = (sqrt((rho+1)^2 + z^2) + sqrt((rho-1)^2 + z^2))/2. It is -1
%   on the crack faces, (2/pi)[1/sqrt(rho^2-1) - asin(1/rho)] ahead of the
%   front in the crack plane, and is independent of Poisson's ratio.
z = abs(z);
l2 = 0.5 * (sqrt((rho + 1).^2 + z.^2) + sqrt((rho - 1).^2 + z.^2));
w = complex(rho.^2 + (z + 1i).^2);
Q = 1 ./ sqrt(w);
on = z == 0 & rho < 1;               % principal branch from z -> 0+
Q(on) = -1i ./ sqrt(1 - rho(on).^2);
G = asin(min(1 ./ l2, 1)) - real(Q);
Gz = imag(Q) + real((z + 1i) .* Q.^3);
s = -(2 / pi) * (G - z .* Gz);
s(on) = -1;
s(z == 0 & rho == 1) = Inf;
end
