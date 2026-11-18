function pop = defectPopulation(type, varargin)
%FLPM.DEFECTPOPULATION  Statistical descriptor of a defect population.
%
%   Instead of reading individual defects, the model is driven by population
%   descriptors (Mater. Des. 266 (2026) 116162, Secs. 3-5). Defect sets are
%   drawn from them by flpm.samplePopulation / flpm.monteCarlo.
%
%   VOLUMETRIC (pores, lack of fusion; XCT)
%     pop = flpm.defectPopulation('volumetric', ...
%             'Size',   makedist('Weibull','a',0.055,'b',5), ... % sqrt(area), mm
%             'Radius', radialData, ...                         % mm from axis
%             'Aspect', aspectData, ...                         % a/c <= 1
%             'Count',  12);                                    % per section
%
%   SURFACE (roughness valleys; profilometry)
%     pop = flpm.defectPopulation('surface', ...
%             'Size', makedist('Exponential','mu',0.01), ...    % Rv - u, mm
%             'Threshold', 0.020, 'SizeOffset', 0.020, ...      % EVT tail
%             'CountPerLength', 0.37, ...                       % per mm
%             'ShieldingFactor', 1.02);                         % Y_R
%
%   Every descriptor accepts a distribution object (makedist/fitdist,
%   including 'Kernel'), a data vector (kernel fit), a scalar (constant) or
%   a makedist cell such as {'GeneralizedPareto','k',-0.29,'sigma',0.0147,
%   'theta',0.050}.
%
%   Name-value pairs
%     'Size'            sqrt(area) (volumetric) or valley depth Rv (surface),
%                       mm. Required.
%     'SizeOffset'      constant added to every sampled size (e.g. the EVT
%                       threshold u when 'Size' describes exceedances)  [0]
%     'Radius'          radial position of the defect centre, mm (volumetric)
%     'Depth'           alternatively, depth of the centre below the free
%                       surface, mm (radius = D/2 - depth)
%     'Aspect'          projected aspect ratio a/c (volumetric)          [1]
%     'Angle'           angular position, rad          [uniform on (-pi, pi)]
%     'Count'           number of defects per cross-section, i.e. the
%                       normalised defect number N (Eq. 5) for volumetric or
%                       V*pi*D (Eq. 6) for surface populations
%     'CountPerLength'  surface only: V, defects per mm of surface line; the
%                       count is round(V*pi*D) for the specimen geometry
%     'CountModel'      'fixed' (published) | 'poisson' (count ~ Poisson with
%                       mean Count)                                  ['fixed']
%     'Threshold'       EVT threshold u (mm): sampled defects with
%                       size <= u are discarded (peaks-over-threshold
%                       reduction, Mater. Des. Sec. 5.4)                [none]
%     'SizeBounds', 'RadiusBounds', 'AspectBounds'
%                       truncation of the descriptors [defaults: [1e-6,Inf],
%                       [0,D/2], (0,1]]
%     'ShieldingFactor' surface only: Y_R, the roughness shielding factor
%                       (Mater. Des. Table 4; see flpm.roughnessShielding)
%                       [geometry SurfaceY]
%     'Name'            label used in results                   [type]
%
%   See also flpm.samplePopulation, flpm.monteCarlo, flpm.roughnessShielding.

type = validatestring(type, {'volumetric', 'surface'});
ip = inputParser;
ip.FunctionName = 'flpm.defectPopulation';
ip.addParameter('Size', []);
ip.addParameter('SizeOffset', 0, @(x) isscalar(x) && x >= 0);
ip.addParameter('Radius', []);
ip.addParameter('Depth', []);
ip.addParameter('Aspect', 1);
ip.addParameter('Angle', makedist('Uniform', 'lower', -pi, 'upper', pi));
ip.addParameter('Count', []);
ip.addParameter('CountPerLength', []);
ip.addParameter('CountModel', 'fixed');
ip.addParameter('Threshold', [], @(x) isempty(x) || (isscalar(x) && x >= 0));
ip.addParameter('SizeBounds', [1e-6 Inf]);
ip.addParameter('RadiusBounds', []);
ip.addParameter('AspectBounds', [0 1]);
ip.addParameter('ShieldingFactor', []);
ip.addParameter('Name', type);
ip.parse(varargin{:});
r = ip.Results;

if isempty(r.Size)
    error('flpm:population:size', 'A ''Size'' descriptor is required.');
end
pop.type = type;
pop.name = char(r.Name);
pop.size = flpm.internal.toDistribution(r.Size, 'Size');
pop.sizeOffset = r.SizeOffset;
pop.sizeBounds = r.SizeBounds;
pop.angle = flpm.internal.toDistribution(r.Angle, 'Angle');
pop.countModel = validatestring(r.CountModel, {'fixed', 'poisson'});
pop.threshold = r.Threshold;

switch type
    case 'volumetric'
        if isempty(r.Radius) == isempty(r.Depth)
            error('flpm:population:location', ...
                'Give exactly one of ''Radius'' or ''Depth'' for a volumetric population.');
        end
        if ~isempty(r.Radius)
            pop.location = 'radius';
            pop.position = flpm.internal.toDistribution(r.Radius, 'Radius');
        else
            pop.location = 'depth';
            pop.position = flpm.internal.toDistribution(r.Depth, 'Depth');
        end
        pop.positionBounds = r.RadiusBounds;
        pop.aspect = flpm.internal.toDistribution(r.Aspect, 'Aspect');
        pop.aspectBounds = r.AspectBounds;
        if isempty(r.Count)
            error('flpm:population:count', 'Give ''Count'' for a volumetric population.');
        end
        pop.count = r.Count;
        pop.countPerLength = [];
        pop.shielding = [];
    case 'surface'
        pop.location = 'surface';
        pop.position = [];
        pop.positionBounds = [];
        pop.aspect = flpm.internal.toDistribution(1, 'Aspect');
        pop.aspectBounds = [0 1];
        if isempty(r.Count) == isempty(r.CountPerLength)
            error('flpm:population:count', ...
                'Give exactly one of ''Count'' or ''CountPerLength'' for a surface population.');
        end
        pop.count = r.Count;
        pop.countPerLength = r.CountPerLength;
        pop.shielding = r.ShieldingFactor;
end
if ~isempty(pop.count) && ~(isscalar(pop.count) && pop.count >= 0)
    error('flpm:population:count', '''Count'' must be a non-negative scalar.');
end
end
