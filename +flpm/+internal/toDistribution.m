function d = toDistribution(spec, name)
%TODISTRIBUTION  Normalise a descriptor specification to a distribution.
%   Accepted specifications
%     prob.ProbabilityDistribution  any object from makedist/fitdist
%                                   (Weibull, GeneralizedPareto, Lognormal,
%                                   Kernel, ...), used as is
%     scalar                        constant value
%     vector of observations        kernel density fit (fitdist 'kernel'),
%                                   as used in Mater. Des. 2026, Sec. 3
%     cell {'Name', p1, v1, ...}    passed to makedist, e.g.
%                                   {'Weibull','a',0.06,'b',4}
%   Returns a struct with fields kind ('dist' | 'const') and value.
if nargin < 2, name = 'descriptor'; end
if isa(spec, 'prob.ProbabilityDistribution')
    d = struct('kind', 'dist', 'value', spec);
elseif isnumeric(spec) && isscalar(spec)
    d = struct('kind', 'const', 'value', double(spec));
elseif isnumeric(spec) && isvector(spec)
    x = double(spec(:));
    x = x(isfinite(x));
    if numel(x) < 2
        error('flpm:descriptor:data', '%s: need at least two observations.', name);
    end
    d = struct('kind', 'dist', 'value', fitdist(x, 'Kernel'));
elseif iscell(spec) && ~isempty(spec) && (ischar(spec{1}) || isstring(spec{1}))
    d = struct('kind', 'dist', 'value', makedist(spec{:}));
else
    error('flpm:descriptor:type', ...
        '%s: expected a distribution object, scalar, data vector or makedist cell.', name);
end
end
