function x = descriptorIcdf(d, u, bounds)
%DESCRIPTORICDF  Inverse CDF of a descriptor truncated to bounds, at u in (0,1).
if strcmp(d.kind, 'const')
    x = repmat(d.value, size(u));
    return
end
lo = bounds(1); hi = bounds(2);
Flo = 0; Fhi = 1;
if isfinite(lo), Flo = cdf(d.value, lo); end
if isfinite(hi), Fhi = cdf(d.value, hi); end
if Fhi - Flo <= 0
    error('flpm:sample:bounds', ...
        'Descriptor has no probability mass inside [%g, %g].', lo, hi);
end
x = icdf(d.value, Flo + u * (Fhi - Flo));
x = min(max(x, lo), hi);              % guard against round-off at the ends
end
