function Y = alegreEmbeddedY(a, aspect, h, Rbar)
%ALEGREEMBEDDEDY  Geometry factor of an embedded elliptical crack in a bar.
%   Y = flpm.internal.alegreEmbeddedY(a, aspect, h, Rbar) returns the larger
%   of the stress-intensity factors at the minor-axis vertex nearer the free
%   surface (Fa1) and at the major-axis vertex (Fc), following the closed-
%   form equations of Alegre, Cuesta & Diaz, Theor. Appl. Fract. Mech. 121
%   (2022) 103438, exactly as transcribed in the published code
%   (closed_form_stress_intensity.m). a and h in the same unit as Rbar.
%   The validity limits of the fit are enforced by clamping
%   0.05 <= (a+h)/R <= 1 and 0.05 <= a/(a+h) <= 0.95 (Fa1 only, as in the
%   published code).
%
%   NOTE (see docs/REVIEW.md): the transcription of M0 contains two lambda^2
%   terms (3.1480*lam^2 ... +3.8745*lam^2); one of them is presumably the
%   lambda^4 term. It is kept as published until checked against the source.
lam = (a + h) ./ Rbar;
lam = min(max(lam, 0.05), 1);

M0 = 3.1480.*lam.^2 - 6.1104.*lam.^3 + 3.8745.*lam.^2 + 1.2555.*lam + 0.6319;
M1 = -1.688.*lam.^4 + 1.2160.*lam.^3 + 1.7455.*lam.^2 - 4.3137.*lam - 2.057;
M2 = -1.6458.*lam.^4 + 4.7320.*lam.^3 - 4.9173.*lam.^2 + 3.1147.*lam - 0.03859;
C1 = M2.*aspect.^2 + M1.*aspect + M0;
B1 = 0.0691.*aspect.^2 - 0.4728.*aspect + 1.0416;
m1 = 1.74.*aspect + (2.0979.*lam.^3 - 2.4414.*lam.^2 - 1.1883.*lam + 3.7175);
r1 = a ./ (a + h);
r1 = min(max(r1, 0.05), 0.95);
Fa1 = B1 ./ (1 - C1 .* r1.^m1);

G0 = 51.3498.*lam.^4 - 79.5268.*lam.^3 + 42.74.*lam.^2 - 4.5145.*lam + 0.2993;
G1 = -81.7186.*lam.^4 + 118.4910.*lam.^3 - 57.1224.*lam.^2 + 2.9721.*lam - 0.1924;
G2 = 32.9938.*lam.^4 - 43.1064.*lam.^3 + 16.8509.*lam.^2 + 1.3102.*lam - 0.0249;
C3 = G2.*aspect.^2 + G1.*aspect + G0;
B3 = 0.5135.*aspect.^3 - 1.4051.*aspect.^2 + 1.3199.*aspect + 0.2105;
m3 = 1.6287.*lam.^3 - 2.8580.*lam.^2 + 2.7763.*lam + 2.3625;
Fc = B3 ./ (1 - C3 .* (a ./ (a + h)).^m3);

% The published code takes max(real(.)) of both; complex values arise only
% if the base of a fractional power is negative, which cannot happen for
% a, h > 0.
Y = max(real(Fa1), real(Fc));
end
