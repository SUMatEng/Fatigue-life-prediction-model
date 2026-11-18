classdef FlpmPhysicsTest < matlab.unittest.TestCase
    %FLPMPHYSICSTEST  Unit tests of the constitutive pieces of the model.
    %   Run with:  results = runtests('tests');

    methods (Test)
        function newmanClosureRatio(tc)
            U = flpm.internal.newmanU(0.1, 2.73);
            tc.verifyEqual(U, 0.782309004496420, 'AbsTol', 1e-12);
            tc.verifyLessThanOrEqual(flpm.internal.newmanU(0.7, 2.73), 1);
        end

        function murakamiInternalAndSurface(tc)
            a = 0.05e-3;  r = a/sqrt(pi);
            S = state([a a a], [1e-3, r, 0]);            % deep, touching, centred
            [Y, A, cls] = flpm.internal.stressIntensity(S);
            tc.verifyEqual(cls(:).', [0 2 2]);
            tc.verifyEqual(Y(:).', [0.5 0.65 0.65]);
            tc.verifyEqual(A(1), a, 'RelTol', 1e-12);
            tc.verifyEqual(A(2), a, 'RelTol', 1e-12);        % whole disc below
            tc.verifyEqual(A(3), a/sqrt(2), 'RelTol', 1e-12); % half disc
        end

        function nearSurfaceRuleAndContinuity(tc)
            a = 0.05e-3;  r = a/sqrt(pi);
            zc = r + r/0.8;                                % proximity limit
            S = state([a a], [zc*(1-1e-9), zc*(1+1e-9)]);
            [Y, A, cls] = flpm.internal.stressIntensity(S);
            tc.verifyEqual(cls(:).', [1 0]);
            tc.verifyEqual(A(1), sqrt(pi*r*(zc + r)/2), 'RelTol', 1e-6);
            tc.verifyEqual(Y(1), 0.65);
            % near-surface envelope equals the disc when it touches the surface
            S = state(a, r*(1+1e-12));
            [~, A] = flpm.internal.stressIntensity(S);
            tc.verifyEqual(A, a, 'RelTol', 1e-9);
        end

        function surfaceAreaDecreasesAsCentreRises(tc)
            a = 0.05e-3; r = a/sqrt(pi);
            z = linspace(r, 0, 20);
            [~, A] = flpm.internal.stressIntensity(state(a*ones(1,20), z));
            tc.verifyTrue(all(diff(A) < 0));
        end

        function noSpecimenDimension(tc)
            % lateral translation of the whole set leaves the life unchanged
            mat = flpm.material('Ti64-AB');
            d1 = flpm.defectSet([0.06 0.05], [0.2 0.25], [0 0.12]);
            d2 = flpm.defectSet([0.06 0.05], [0.2 0.25], [0 0.12] + 1234.5);
            r1 = flpm.fatigueLife(300, d1, mat);
            r2 = flpm.fatigueLife(300, d2, mat);
            tc.verifyEqual(r2.N, r1.N, 'RelTol', 1e-12);
        end

        function roughnessShieldingFactor(tc)
            [YR, sa] = flpm.roughnessShielding(0.010, 0.100);
            tc.verifyEqual(YR, 1.032170528546519, 'AbsTol', 1e-12);
            tc.verifyEqual(sa, 0.1*(2.97*0.1 - 3.51*0.01 - 9.74*0.001), 'AbsTol', 1e-15);
            [~, sa2] = flpm.roughnessShielding(0.05, 0.1);
            tc.verifyEqual(sa2, 0.038, 'AbsTol', 1e-15);
            tc.verifySize(flpm.roughnessShielding([0.01 0.02 0.05], 0.1), [1 3]);
        end

        function coalescenceMergesByArea(tc)
            a1 = 0.06e-3; a2 = 0.04e-3;
            r = (a1 + a2)/sqrt(pi);
            S = state([a1 a2], [0.5e-3 0.5e-3], [0 0.9*r]);
            S.ext = [2e-6; 5e-6];
            S.dK = [2; 1];
            [M, merged] = flpm.internal.coalesce(S);
            tc.verifyTrue(merged);
            tc.verifyNumElements(M.a, 1);
            tc.verifyEqual(M.a, hypot(a1, a2), 'RelTol', 1e-12);
            tc.verifyEqual(M.s, a2^2*0.9*r/(a1^2 + a2^2), 'RelTol', 1e-12);
            tc.verifyEqual(M.ext, 5e-6);
            tc.verifyEqual(M.id, 1);
        end

        function coalescenceIsTransitive(tc)
            a = 0.05e-3; d = 1.9*a/sqrt(pi);
            S = state([a a a], [1e-3 1e-3 1e-3], [0 d 2*d]);
            S.dK = [1; 3; 2];
            M = flpm.internal.coalesce(S);
            tc.verifyNumElements(M.a, 1);
            tc.verifyEqual(M.id, 2);
            tc.verifyEqual(M.a, sqrt(3)*a, 'RelTol', 1e-12);
        end

        function noCoalescenceWhenApart(tc)
            a = 0.05e-3;
            S = state([a a], [1e-3 1e-3], [0 3*a]);
            [M, merged] = flpm.internal.coalesce(S);
            tc.verifyFalse(merged);
            tc.verifyNumElements(M.a, 2);
        end
    end
end

function S = state(a, z, s)
% Solver state (SI) from vectors of sqrt(area), depth and position.
if nargin < 3, s = zeros(size(a)); end
n = numel(a);
S = struct('a', a(:), 'z', z(:), 's', s(:), 'ext', zeros(n,1), ...
    'id', (1:n)', 'dK', zeros(n,1));
end
