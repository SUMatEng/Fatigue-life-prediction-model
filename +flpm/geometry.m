function geom = geometry(varargin)
%FLPM.GEOMETRY  Specimen geometry and loading for the fatigue-strength model.
%
%   geom = flpm.geometry('Diameter', 7.3, 'LoadRatio', 0.1)
%
%   Name-value pairs (defaults in brackets):
%     'Diameter'      round-bar diameter D in mm                      [7.3]
%     'LoadRatio'     R = Smin/Smax                                    [0.1]
%     'Eccentricity'  load-line eccentricity factor e used in the
%                     net-section stress of a part-through crack       [0.13]
%     'SurfaceY'      lower bound on Y for a defect that breaks the
%                     free surface (roughness shielding factor; the
%                     2023-2025 code used 1.0335 for all such defects) [1.0335]
%
%   Geometry is in mm (as XCT/profilometry data are); the solver works
%   internally in SI units.

ip = inputParser;
ip.FunctionName = 'flpm.geometry';
ip.addParameter('Diameter', 7.3, @(x) isscalar(x) && x > 0);
ip.addParameter('LoadRatio', 0.1, @(x) isscalar(x) && x < 1);
ip.addParameter('Eccentricity', 0.13, @(x) isscalar(x) && x >= 0);
ip.addParameter('SurfaceY', 1.0335, @(x) isscalar(x) && x > 0);
ip.parse(varargin{:});
r = ip.Results;
geom = struct('D', r.Diameter, 'R', r.LoadRatio, ...
    'e', r.Eccentricity, 'surfaceY', r.SurfaceY);
end
