function [sqrtArea, area] = projectedArea(varargin)
%FLPM.PROJECTEDAREA  Murakami sqrt(area) of 3D defects on the plane normal to the load.
%
%   s = flpm.projectedArea('Projections', [Ayz Axz Axy], 'Normal', n)
%   s = flpm.projectedArea('Ellipsoid', [a b c], 'Axes', Q, 'Normal', n)
%   [s, A] = ...                 s = sqrt(A), in the length unit of the input
%
%   Murakami's sqrt(area) is the area of the defect projected onto the
%   plane normal to the maximum principal stress. n is that direction (a
%   3-vector, or one row per defect; for uniaxial loading, the load axis).
%
%   'Projections'  N x 3 projected areas of each defect on the yz, xz and xy
%                  planes (the XCT exports ProjectedArea_yz/xz/xy). The
%                  defect is taken as the ellipsoid, aligned with x, y, z,
%                  with these three projections; for it the projection on
%                  the plane normal to n is exactly
%                      A(n) = sqrt((n_x Ayz)^2 + (n_y Axz)^2 + (n_z Axy)^2).
%                  For n along an axis this returns that projection
%                  unchanged, so it is exact for any shape in that case.
%   'Ellipsoid'    N x 3 semi-axes [a b c] of the equivalent ellipsoid of
%                  each defect (e.g. from its inertia tensor), with
%   'Axes'         3 x 3 x N (or 3 x 3) principal directions as columns
%                  (default: aligned with x, y, z). For the ellipsoid
%                  x' M x <= 1, M = Q diag(a^-2, b^-2, c^-2) Q',
%                      A(n) = pi a b c sqrt(n' M n).
%
%   For tilted or irregular defects away from the axes, the ellipsoid is
%   an approximation of the true silhouette (convex hull of the defect in
%   the plane), which is what the XCT software should export when the load
%   direction is not a scan axis.
%
%   See also flpm.defectSet.

ip = inputParser;
ip.FunctionName = 'flpm.projectedArea';
ip.addParameter('Projections', []);
ip.addParameter('Ellipsoid', []);
ip.addParameter('Axes', eye(3));
ip.addParameter('Normal', [0 0 1]);
ip.parse(varargin{:});
o = ip.Results;
if isempty(o.Projections) == isempty(o.Ellipsoid)
    error('flpm:projectedArea:input', 'Give either ''Projections'' or ''Ellipsoid''.');
end
n = o.Normal;
if size(n, 2) ~= 3, n = reshape(n, 1, 3); end
n = n ./ vecnorm(n, 2, 2);
if ~isempty(o.Projections)
    P = o.Projections;
    if size(P, 2) ~= 3 || any(P(:) < 0)
        error('flpm:projectedArea:input', '''Projections'' must be N x 3 and non-negative.');
    end
    area = sqrt((n(:,1) .* P(:,1)).^2 + (n(:,2) .* P(:,2)).^2 + (n(:,3) .* P(:,3)).^2);
else
    E = o.Ellipsoid;
    if size(E, 2) ~= 3 || any(E(:) <= 0)
        error('flpm:projectedArea:input', '''Ellipsoid'' must be N x 3 and positive.');
    end
    N = size(E, 1);
    Q = o.Axes;
    if size(Q, 3) == 1, Q = repmat(Q, 1, 1, N); end
    if size(n, 1) == 1, n = repmat(n, N, 1); end
    area = zeros(N, 1);
    for k = 1:N
        u = Q(:,:,k)' * n(k,:)';           % n in the principal frame
        area(k) = pi * prod(E(k,:)) * sqrt(sum((u(:)' ./ E(k,:)).^2));
    end
end
sqrtArea = sqrt(area);
end
