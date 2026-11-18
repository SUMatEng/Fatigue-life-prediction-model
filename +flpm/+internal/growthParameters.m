function [C1, dKe, dKl, noClosure] = growthParameters(S, mat)
%GROWTHPARAMETERS  NASGRO C1 and thresholds of every crack (environment).
%   Cracks whose disc does not reach the free surface (z > r) are sealed
%   from the air; they take the constants in mat.sealed where given.
%   mat.sealed.closure = false removes closure for them (U = 1).
n = numel(S.a);
C1 = mat.C1 * ones(n, 1);
dKe = mat.dKth_eff * ones(n, 1);
dKl = mat.dKth_lc * ones(n, 1);
sealed = S.z > S.a / sqrt(pi);
noClosure = false(n, 1);
if isempty(mat.sealed) || ~any(sealed), return, end
e = mat.sealed;
if isfield(e, 'C1'), C1(sealed) = e.C1; end
if isfield(e, 'dKth_eff'), dKe(sealed) = e.dKth_eff; end
if isfield(e, 'dKth_lc'), dKl(sealed) = e.dKth_lc; end
if isfield(e, 'closure') && ~e.closure, noClosure = sealed; end
end
