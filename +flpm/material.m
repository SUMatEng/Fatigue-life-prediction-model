function mat = material(name, varargin)
%FLPM.MATERIAL  NASGRO / short-crack parameters for PBF-LB alloys.
%
%   mat = flpm.material(name) returns the parameters of a predefined
%   material. Name-value pairs override fields, e.g.
%       mat = flpm.material('Ti64-AB', 'KIC', 55);
%
%   Predefined materials (Macallister, Deillon & Becker, Mater. Des. 266
%   (2026) 116162, Table 2; Macallister & Becker, Acta Mater. 237 (2022)
%   118189, Table 1):
%       'Ti64-AB'      PBF-LB Ti-6Al-4V, as-built / as-fabricated
%       'Ti64-AN'      PBF-LB Ti-6Al-4V, annealed
%       'AlSi10Mg-AB'  PBF-LB Al-10Si-Mg, as-built
%
%   Fields (MPa, MPa*sqrt(m), m/cycle):
%       C1, p, m      NASGRO constants (Zhang/Maierhofer form)
%       alpha         constraint factor in Newman's closure function
%       dKth_lc       long-crack threshold at the load ratio of the data
%       dKth_eff      intrinsic threshold
%       dSigma_e      defect-free fatigue strength (stress range)
%       KIC           fracture toughness
%       sigma0        cyclic yield strength (plastic-zone size)
%       dKth_cov      coefficient of variation of dKth_lc between specimens
%                     (lognormal, median dKth_lc); 0 = deterministic    [0]
%                     It also scatters the closure length l, which
%                     follows from dKth_lc through the El Haddad length.
%
%   See also flpm.options.

switch lower(strrep(name, '_', '-'))
    case {'ti64-ab', 'ti64-af', 'ti-6al-4v'}
        mat = base('Ti64-AB', 1.61e-10, 0.35, 2.77, 2.73, 2.70, 1.80, 450, 50.0, 1100);
    case {'ti64-an'}
        mat = base('Ti64-AN', 1.87e-11, 0.35, 2.94, 1.11, 5.6, 2.7, 487, 50.0, 900);
    case {'alsi10mg-ab', 'alsi10mg', 'al-10si-mg'}
        mat = base('AlSi10Mg-AB', 12.7e-9, 0.48, 1.30, 3.53, 1.26, 1.03, 182, 37.4, 304);
    otherwise
        error('flpm:material:unknown', 'Unknown material "%s".', name);
end
if mod(numel(varargin), 2) ~= 0
    error('flpm:material:args', 'Overrides must be name-value pairs.');
end
for k = 1:2:numel(varargin)
    f = varargin{k};
    if ~isfield(mat, f)
        error('flpm:material:field', 'Unknown material field "%s".', f);
    end
    mat.(f) = varargin{k+1};
end
end

function mat = base(name, C1, p, m, alpha, dKlc, dKeff, dSe, KIC, s0)
mat = struct('name', name, 'C1', C1, 'p', p, 'm', m, 'alpha', alpha, ...
    'dKth_lc', dKlc, 'dKth_eff', dKeff, 'dSigma_e', dSe, 'KIC', KIC, ...
    'sigma0', s0, 'dKth_cov', 0);
end
