function m = thresholdDraw(mat, seed, k)
%THRESHOLDDRAW  Material of simulation k with a random dKth_lc.
%   Uses its own stream (mlfg6331_64, substream k), so the defect sampling
%   of flpm.samplePopulation is unaffected.
m = mat;
if ~isfield(mat, 'dKth_cov') || mat.dKth_cov <= 0, return, end
rs = RandStream('mlfg6331_64', 'Seed', seed);
rs.Substream = k;
s = sqrt(log(1 + mat.dKth_cov^2));
m.dKth_lc = max(mat.dKth_lc * exp(s * randn(rs)), mat.dKth_eff * (1 + 1e-6));
m.dKth_cov = 0;
end
