function opts = modelOptions(preset, varargin)
%FLPM.MODELOPTIONS  Physics switches and numerical settings for the solver.
%
%   opts = flpm.modelOptions()            % 'paper' preset (default)
%   opts = flpm.modelOptions('legacy')    % reproduces the 2023-2025 code
%   opts = flpm.modelOptions('paper', 'TargetLife', 2e6, 'Eta', 0.01)
%
%   PRESETS
%   'paper'   Equations as written in Macallister & Becker (Acta Mater. 2022)
%             and Macallister, Deillon & Becker (Mater. Des. 2026), solved
%             with adaptive cycle stepping and a bracketed root search.
%             Differs from the legacy code only where the code disagrees
%             with the papers (see docs/REVIEW.md):
%               - plastic zone rp = (1/6pi)(Kmax/sigma0)^2 (Acta Eq. 12);
%                 the code used (pi/6)(Kmax/sigma0)^2, i.e. pi^2 larger.
%               - aspect-ratio update applied once per increment (the code
%                 applied the increment twice) and bounded to
%                 [0.2, 1] (Mater. Des. Sec. 2, "initial minimum of 0.2").
%   'legacy'  Reproduces the numerics of Nasgro.m / multi_coalition_version3.m
%             (N. Macallister, Feb 2025), including fixed 100-cycle steps and
%             the fminsearch objective, but without the crashes. Use it
%             for regression against the published results.
%
%   OPTIONS (name-value, override any preset value)
%   Physics
%     PlasticZoneFactor   c in rp = c (Kmax/sigma0)^2 [1/(6*pi) | pi/6]
%     AspectUpdate        'paper' | 'legacy'
%     AspectBounds        [lo hi] bounds on a/c after each update
%     CreepRate           da/dN (m/cycle) assigned to a defect whose
%                         effective dK exceeds dK_th,eff but whose dK is below
%                         the current cyclic R-curve threshold. 1e-12 in the
%                         published model; 0 gives true crack arrest.  [1e-12]
%     InitiationStress    stress used in the growth (initiation) check:
%                         'nominal' (published code) | 'net' (net-section
%                         stress, as described in Acta Sec. 2.2)    ['nominal']
%     InteractionRadius   circle radius used for coalescence and the
%                         surface proximity rule: 'sqrtArea' (published code,
%                         r = sqrt(area)) | 'equivalentCircle'
%                         (r = sqrt(area/pi))                    ['sqrtArea']
%     NetSectionInDeltaK  use the growing net-section stress inside dK
%                         (published model) [true]
%     NetSectionKCap      net-section stress is only updated while
%                         max(U*dK) is below this value (MPa*sqrt(m))  [20]
%     InternalYFactor     factor on Y for internal defects       [0.5/0.65]
%     NearSurfaceRatio    a/h above which a defect is near-surface  [0.8]
%   Numerics
%     TargetLife          cycles defining the fatigue strength     [1e7]
%     MaxLife             cycles beyond which a run is a run-out
%                         ['paper': 1.02*TargetLife, 'legacy': TargetLife+2e5]
%     Stepping            'adaptive' | 'fixed'
%     FixedStep           cycles per increment when Stepping='fixed' [100]
%     Eta                 adaptive step: max crack advance per increment as
%                         a fraction of min(a, closure length l, gaps) [0.005]
%     MaxStep, MinStep    bounds on the adaptive increment (cycles) [1e4, 1]
%     Solver              'bracket' (bisection/secant on log N) |
%                         'fminsearch' (legacy objective)
%     StressTol           absolute stress tolerance, MPa (bracket)  [0.05]
%     InitialStress       start value for fminsearch, MPa           [400]
%     RecordHistory       store crack-growth history in fatigueLife [false]

if nargin < 1 || isempty(preset)
    preset = 'paper';
end

switch lower(preset)
    case 'paper'
        opts = struct( ...
            'preset', 'paper', ...
            'PlasticZoneFactor', 1/(6*pi), ...
            'AspectUpdate', 'paper', ...
            'AspectBounds', [0.2 1], ...
            'CreepRate', 1e-12, ...
            'InitiationStress', 'nominal', ...
            'InteractionRadius', 'sqrtArea', ...
            'NetSectionInDeltaK', true, ...
            'NetSectionKCap', 20, ...
            'InternalYFactor', 0.5/0.65, ...
            'NearSurfaceRatio', 0.8, ...
            'TargetLife', 1e7, ...
            'MaxLife', [], ...
            'Stepping', 'adaptive', ...
            'FixedStep', 100, ...
            'Eta', 0.005, ...
            'MaxStep', 1e4, ...
            'MinStep', 1, ...
            'Solver', 'bracket', ...
            'StressTol', 0.05, ...
            'InitialStress', 400, ...
            'RecordHistory', false);
    case 'legacy'
        opts = flpm.modelOptions('paper');
        opts.preset = 'legacy';
        opts.PlasticZoneFactor = pi/6;
        opts.AspectUpdate = 'legacy';
        opts.AspectBounds = [0 Inf];
        opts.Stepping = 'fixed';
        opts.Solver = 'fminsearch';
        opts.MaxLife = [];      % resolved below to TargetLife + 2e5
    otherwise
        error('flpm:modelOptions:preset', 'Unknown preset "%s".', preset);
end

if mod(numel(varargin), 2) ~= 0
    error('flpm:modelOptions:args', 'Options must be name-value pairs.');
end
for k = 1:2:numel(varargin)
    f = varargin{k};
    match = strcmpi(fieldnames(opts), f);
    if ~any(match)
        error('flpm:modelOptions:field', 'Unknown option "%s".', f);
    end
    fn = fieldnames(opts);
    opts.(fn{match}) = varargin{k+1};
end

if isempty(opts.MaxLife)
    if strcmp(opts.preset, 'legacy')
        opts.MaxLife = opts.TargetLife + 2e5;
    else
        opts.MaxLife = 1.02*opts.TargetLife;
    end
end
end
