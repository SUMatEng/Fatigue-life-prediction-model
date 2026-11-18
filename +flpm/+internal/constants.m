function k = constants()
%CONSTANTS  Fixed constants of the model (Murakami sqrt(area) framework).
%   Yinternal   geometry factor of an internal defect, K = 0.5 s sqrt(pi sqrt(area))
%   Ysurface    geometry factor of a surface or near-surface defect (0.65)
%   proximity   Murakami's rule: an internal defect whose ligament to the
%               free surface is h <= r/0.8 (r = equivalent-circle radius)
%               is treated as a surface defect enclosing the ligament
%   rpFactor    cyclic plastic-zone size rp = rpFactor (Kmax/sigma0)^2,
%               Acta Mater. 237 (2022) Eq. 12
k.Yinternal = 0.5;
k.Ysurface  = 0.65;
k.proximity = 0.8;
k.rpFactor  = 1 / (6*pi);
end
