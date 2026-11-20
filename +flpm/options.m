function opts = options(varargin)
%FLPM.OPTIONS  Loading and numerical settings.
%
%   opts = flpm.options('LoadRatio', 0.1, 'TargetLife', 1e7)
%
%   These are the only settings; the physics has a single form (see
%   docs/MODEL.md).
%
%   'LoadRatio'      R = sigma_min/sigma_max                         [0.1]
%   'TargetLife'     N at which the fatigue strength is defined      [1e7]
%   'Eta'            adaptive increment: largest relative advance of
%                    any crack per increment                        [0.005]
%   'MaxStep', 'MinStep'  bounds on the increment (cycles)       [1e4, 1]
%   'StressTol'      stress tolerance of the strength solver (MPa)  [0.05]
%   'RecordHistory'  keep N, a_max, number of defects and crack
%                    snapshots (for flpm.plotProgression)      [false]
%   'ResidualStress' residual stress normal to the crack plane (MPa,
%                    tensile positive) as a function of depth below the
%                    free surface (mm), e.g. @(z) -400*exp(-z/0.1);
%                    [] = none. It shifts the local mean stress: at a
%                    defect centred at depth z,
%                      Kmax = (sigma_max + sigma_res(z)) Y sqrt(pi A),
%                      R_loc = (sigma_min + sigma_res)/(sigma_max + sigma_res),
%                    with the range unchanged and U = U(R_loc). A crack
%                    with sigma_max + sigma_res <= 0 is closed through the
%                    whole cycle and does not grow. Evaluated at the
%                    centre of each defect (no gradient over the crack,
%                    no relaxation).                                  [[]]
ip = inputParser;
ip.FunctionName = 'flpm.options';
ip.addParameter('LoadRatio', 0.1, @(x) isscalar(x) && x < 1);
ip.addParameter('TargetLife', 1e7, @(x) isscalar(x) && x > 0);
ip.addParameter('Eta', 0.005, @(x) isscalar(x) && x > 0 && x < 1);
ip.addParameter('MaxStep', 1e4, @(x) isscalar(x) && x >= 1);
ip.addParameter('MinStep', 1, @(x) isscalar(x) && x > 0);
ip.addParameter('StressTol', 0.05, @(x) isscalar(x) && x > 0);
ip.addParameter('RecordHistory', false, @(x) islogical(x) || isnumeric(x));
ip.addParameter('ResidualStress', [], @(x) isempty(x) || isa(x, 'function_handle'));
ip.parse(varargin{:});
opts = ip.Results;
opts.RecordHistory = logical(opts.RecordHistory);
% A run is a run-out once it exceeds the target life by 2 %, which the
% strength solver then brackets.
opts.MaxLife = 1.02 * opts.TargetLife;
end
