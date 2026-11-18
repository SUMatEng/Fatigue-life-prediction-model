function c = toCopula(spec)
%TOCOPULA  Normalise a size-depth copula specification.
%   {'Gaussian', rho} | {'t', rho, nu} | struct with family, rho(, nu) |
%   output of flpm.testDependence (its field 'copula').
if isempty(spec), c = []; return, end
if isstruct(spec) && isfield(spec, 'copula'), spec = spec.copula; end
if isstruct(spec)
    c = spec;
elseif iscell(spec)
    c.family = lower(spec{1});
    c.rho = spec{2};
    c.nu = Inf;
    if numel(spec) > 2, c.nu = spec{3}; end
else
    error('flpm:population:copula', 'Unrecognised ''Copula'' specification.');
end
if isempty(c), return, end
c.family = lower(c.family);
if ~any(strcmp(c.family, {'gaussian', 't'})) || ~(abs(c.rho) < 1)
    error('flpm:population:copula', 'Copula must be Gaussian or t with |rho| < 1.');
end
if strcmp(c.family, 't') && ~(isfield(c, 'nu') && c.nu > 0 && isfinite(c.nu))
    error('flpm:population:copula', 'A t copula needs nu > 0.');
end
end
