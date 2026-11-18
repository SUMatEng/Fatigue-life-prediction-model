function U = newmanU(R, alpha, smaxSigma0)
%NEWMANU  Newman (1981/1984) plasticity-induced closure ratio U = dKeff/dK.
%   U = flpm.internal.newmanU(R, alpha, smaxSigma0)
%   alpha is the constraint factor; smaxSigma0 = Smax/sigma0. If smaxSigma0
%   is empty the legacy rule is used (0.2 for R > 0, 0.3 for R <= 0).
%   For R > 0 the cubic form is used; for R <= 0 the linear form, exactly as
%   in the published code.
% Vectorised over R (one value per defect when residual stress shifts the
% local load ratio).
if nargin < 3 || isempty(smaxSigma0)
    smaxSigma0 = 0.2 * (R > 0) + 0.3 * (R <= 0);
end
% Coefficients exactly as implemented in the published code. Note that
% Newman's (1984) closure function uses A0 = (0.825 - 0.34*alpha +
% 0.05*alpha^2)*cos(pi/2*Smax/sigma0)^(1/alpha); the source of the modified
% coefficients 0.8613/0.3387/0.0465 is not documented (see docs/REVIEW.md).
A0 = (0.8613 - 0.3387*alpha + 0.0465*alpha^2) .* cos(pi/2*smaxSigma0).^(1/alpha);
A1 = (1.047 - 0.233*alpha) .* smaxSigma0;
A3 = 2*A0 + A1 - 1;
A2 = 1 - A0 - A1 - A3;
fPos = max(R, A0 + A1.*R + A2.*R.^2 + A3.*R.^3);
fNeg = max(R, A0 + A1.*R);
f = fNeg;
f(R > 0) = fPos(R > 0);
U = (1 - f) ./ (1 - R);
end
