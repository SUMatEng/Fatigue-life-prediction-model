function out = feCrackK(cracks, varargin)
%FECRACKK  Mode-I K along the fronts of embedded penny cracks by 3D FE.
%
%   out = feCrackK(cracks)
%   cracks: N x 4 [x y z r] (any unit), every crack normal to the load (z).
%
%   Each crack is a thin oblate spheroidal cavity (half-thickness
%   t = 'Thickness' * r) in a block of edge 'Block' * (extent of the
%   cluster), under unit remote tension along z (PDE Toolbox, linear
%   elastic, linear tetrahedra). K at each of 'NumAngles' front points
%   follows from the crack-opening displacement by displacement
%   extrapolation:
%       K_app(d) = E' delta(d) / (8 sqrt(d / 2 pi)),  E' = E/(1 - nu^2)
%   fitted linearly in d over d = r - rho in 'FitRange' * r and taken at
%   d = 0. The ratio to the same crack alone on the same mesh settings
%   (out.gamma) is the interaction factor; out.Kref is the single-crack
%   value, to compare with Sneddon's 2 sqrt(r/pi).
%
%   Name-value pairs: 'Thickness' [0.02], 'MeshFace' [0.03] (element size
%   on the cavities, relative to r), 'NumAngles' [24], 'FitRange'
%   [0.04 0.2], 'Block' [5], 'Reference' [true].
ip = inputParser;
ip.addParameter('Thickness', 0.02);
ip.addParameter('MeshFace', 0.03);
ip.addParameter('NumAngles', 24);
ip.addParameter('FitRange', [0.04 0.2]);
ip.addParameter('Block', 5);
ip.addParameter('Reference', true);
ip.parse(varargin{:});
o = ip.Results;
E = 114e3; nu = 0.34;
c = cracks;
rmin = min(c(:, 4));
out.K = solveSet(c, o, E, nu);
if o.Reference
    out.Kref = zeros(size(out.K));
    for i = 1:size(c, 1)
        ci = c(i, :);
        Ki = solveSet(ci, o, E, nu);
        out.Kref(i, :) = Ki;
    end
    out.gamma = out.K ./ out.Kref;
end
out.cracks = c;
out.theta = (0:o.NumAngles-1) / o.NumAngles * 2 * pi;
out.options = o;
out.Ksneddon = 2 * sqrt(c(:, 4) / pi);
end

function K = solveSet(c, o, E, nu)
n = size(c, 1);
ext = max(max(abs(c(:, 1:3) - mean(c(:, 1:3), 1)) + c(:, 4), [], 1));
L = 2 * o.Block * max(ext, max(c(:, 4)));
ctr = mean(c(:, 1:3), 1);
B = fegeometry(multicuboid(L, L, L));
B = translate(B, ctr - [0 0 L/2]);
for i = 1:n
    G0 = cavity(c(i, :), o.Thickness * c(i, 4));
    B = subtract(B, G0);
end
G = B;
V = G.Vertices;
zt = max(V(:, 3)); zb = min(V(:, 3));
ftop = nearestFace(G, [ctr(1:2) zt]);
fbot = nearestFace(G, [ctr(1:2) zb]);
fc = arrayfun(@(i) nearestFace(G, c(i, 1:3) + [0 0 o.Thickness * c(i, 4)]), 1:n);
bot = find(abs(V(:, 3) - zb) < 1e-9 * L);
[~, k1] = min(V(bot, 1) + V(bot, 2));
[~, k2] = max(V(bot, 1) - V(bot, 2));
model = femodel(AnalysisType = 'structuralStatic', Geometry = G);
model.MaterialProperties = materialProperties(YoungsModulus = E, PoissonsRatio = nu);
model.FaceLoad(ftop) = faceLoad(SurfaceTraction = [0 0 1]);
model.FaceBC(fbot) = faceBC(ZDisplacement = 0);
model.VertexBC(bot(k1)) = vertexBC(XDisplacement = 0, YDisplacement = 0);
model.VertexBC(bot(k2)) = vertexBC(YDisplacement = 0);
hf = arrayfun(@(i) o.MeshFace * c(i, 4), 1:n, 'UniformOutput', false);
hfc = [num2cell(fc); hf];
% retry with other size limits: close crack pairs can fail to mesh
tries = [6 0.15 1.4; 10 0.08 1.4; 8 0.3 1.3; 12 0.05 1.5];
for k = 1:size(tries, 1)
    try
        model = generateMesh(model, 'Hmax', L / tries(k, 1), ...
            'Hmin', tries(k, 2) * o.MeshFace * min(c(:, 4)), ...
            'Hface', hfc(:)', 'Hgrad', tries(k, 3));
        break
    catch err
        if k == size(tries, 1), rethrow(err), end
    end
end
R = solve(model);
P = model.Geometry.Mesh.Nodes;
uz = R.Displacement.uz(:)';
th = (0:o.NumAngles-1) / o.NumAngles * 2 * pi;
K = zeros(n, numel(th));
Ep = E / (1 - nu^2);
for i = 1:n
    r = c(i, 4); t = o.Thickness * r;
    dx = P(1, :) - c(i, 1); dy = P(2, :) - c(i, 2); dz = P(3, :) - c(i, 3);
    on = abs(hypot(dx, dy).^2 / r^2 + dz.^2 / t^2 - 1) < 2e-3;
    up = on & dz > 0; dn = on & dz < 0;
    Fu = scatteredInterpolant(dx(up)', dy(up)', uz(up)', 'natural', 'nearest');
    Fd = scatteredInterpolant(dx(dn)', dy(dn)', uz(dn)', 'natural', 'nearest');
    d = linspace(o.FitRange(1), o.FitRange(2), 12) * r;
    for k = 1:numel(th)
        rho = r - d;
        delta = Fu(rho * cos(th(k)), rho * sin(th(k))) - Fd(rho * cos(th(k)), rho * sin(th(k)));
        Ka = Ep * delta ./ (8 * sqrt(d / (2 * pi)));
        p = polyfit(d, Ka, 1);
        K(i, k) = p(2);
    end
end
end

function G = cavity(ci, t)
[X, Y, Z] = ellipsoid(ci(1), ci(2), ci(3), ci(4), ci(4), t, 120);
F = surf2patch(X, Y, Z, 'triangles');
[V, ~, ic] = uniquetol(F.vertices, 1e-9, 'ByRows', true);
Fc = ic(F.faces);
Fc(Fc(:,1) == Fc(:,2) | Fc(:,2) == Fc(:,3) | Fc(:,1) == Fc(:,3), :) = [];
G = fegeometry(triangulation(Fc, V));
end
