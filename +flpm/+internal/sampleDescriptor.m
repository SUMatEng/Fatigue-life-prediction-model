function x = sampleDescriptor(d, n, bounds, rs)
%SAMPLEDESCRIPTOR  Draw n values from a descriptor, truncated to bounds.
%   Truncation is exact inverse-CDF sampling on [F(lo), F(hi)], so no
%   probability mass is piled up at the bounds (the 2023-2025 scripts
%   clipped out-of-range kernel samples to the bound instead).
if strcmp(d.kind, 'const')
    x = repmat(d.value, n, 1);
    return
end
x = flpm.internal.descriptorIcdf(d, rand(rs, n, 1), bounds);
end
