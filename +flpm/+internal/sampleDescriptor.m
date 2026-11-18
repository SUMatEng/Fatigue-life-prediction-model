function x = sampleDescriptor(d, n, bounds, rs)
%SAMPLEDESCRIPTOR  Draw n values from a descriptor, truncated to bounds.
%   Truncation is exact inverse-CDF sampling on [F(lo), F(hi)], so no
%   probability mass is piled up at the bounds (the 2023-2025 scripts
%   clipped out-of-range kernel samples to the bound instead).
if strcmp(d.kind, 'const')
    x = repmat(d.value, n, 1);
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
u = rand(rs, n, 1);
x = icdf(d.value, Flo + u * (Fhi - Flo));
x = min(max(x, lo), hi);              % guard against round-off at the ends
end
