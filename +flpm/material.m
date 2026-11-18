function mat = material(name, varargin)
%FLPM.MATERIAL  NASGRO / short-crack parameters for PBF-LB alloys.
%
%   mat = flpm.material(name) returns a struct of material parameters for
%   one of the predefined materials. Name-value pairs override fields, e.g.
%       mat = flpm.material('Ti64-AB', 'dKth_lc', 2.4);
%
%   Predefined materials (Macallister, Deillon & Becker, Mater. Des. 266
%   (2026) 116162, Table 2; Macallister & Becker, Acta Mater. 237 (2022)
%   118189, Table 1):
%       'Ti64-AB'      PBF-LB Ti-6Al-4V, as-built / as-fabricated (AF)
%       'Ti64-AN'      PBF-LB Ti-6Al-4V, annealed (Acta Mater. 2022, Table 1)
%       'AlSi10Mg-AB'  PBF-LB Al-10Si-Mg, as-built
%
%   Fields (units: MPa, MPa*sqrt(m), m/cycle):
%       C1, p, m, alpha     NASGRO (Zhang/Maierhofer form) constants
%       dKth_lc             long-crack threshold dK_th,lc at the load ratio R
%       dKth_eff            intrinsic (effective) threshold dK_th,eff
%       dSigma_e            defect-free fatigue strength (range) for El Haddad
%       KIC                 plane-strain fracture toughness
%       sigma0, sigmaU      yield and ultimate tensile strength
%       flowStress          (sigma0+sigmaU)/2 unless overridden (Acta Eq. 7)
%       Yw                  geometry factor used to express El Haddad a0
%                           in sqrt(area) terms (Murakami surface value 0.65)
%       newmanSmaxSigma0    Smax/sigma0 constraint ratio in Newman's U.
%                           [] = legacy rule: 0.2 for R>0, 0.3 for R<=0.
%
%   See also flpm.geometry, flpm.modelOptions.

switch lower(strrep(name, '_', '-'))
    case {'ti64-ab', 'ti64-af', 'ti-6al-4v'}
        mat = base('Ti64-AB', 1.61e-10, 0.35, 2.77, 2.73, 2.70, 1.80, 450, 50.0, 1100, 1200);
    case {'ti64-an'}
        mat = base('Ti64-AN', 1.87e-11, 0.35, 2.94, 1.11, 5.6, 2.7, 487, 50.0, 900, 950);
    case {'alsi10mg-ab', 'alsi10mg', 'al-10si-mg'}
        mat = base('AlSi10Mg-AB', 12.7e-9, 0.48, 1.30, 3.53, 1.26, 1.03, 182, 37.4, 304, 401);
    otherwise
        error('flpm:material:unknown', 'Unknown material "%s".', name);
end

if mod(numel(varargin), 2) ~= 0
    error('flpm:material:args', 'Overrides must be name-value pairs.');
end
flowGiven = false;
for k = 1:2:numel(varargin)
    f = varargin{k};
    if ~isfield(mat, f)
        error('flpm:material:field', 'Unknown material field "%s".', f);
    end
    mat.(f) = varargin{k+1};
    flowGiven = flowGiven || strcmp(f, 'flowStress');
end
if ~flowGiven
    mat.flowStress = 0.5*(mat.sigma0 + mat.sigmaU);
end
end

function mat = base(name, C1, p, m, alpha, dKlc, dKeff, dSe, KIC, s0, su)
mat = struct( ...
    'name', name, ...
    'C1', C1, 'p', p, 'm', m, 'alpha', alpha, ...
    'dKth_lc', dKlc, 'dKth_eff', dKeff, ...
    'dSigma_e', dSe, 'KIC', KIC, ...
    'sigma0', s0, 'sigmaU', su, 'flowStress', 0.5*(s0 + su), ...
    'Yw', 0.65, 'newmanSmaxSigma0', []);
end
