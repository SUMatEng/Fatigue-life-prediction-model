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

        function linkUpAndRecharacterisation(tc)
            % overlapping coplanar flaws link up; replaced by the ellipse
            % inscribed in the containing rectangle (BS 7910)
            a1 = 0.06e-3; a2 = 0.04e-3; r1 = a1/sqrt(pi); r2 = a2/sqrt(pi);
            S = state([a1 a2], [0.5e-3 0.5e-3], [0, r1 + r2 - 0.1*r2]);
            S.ext = [2e-6; 5e-6];  S.dK = [2; 1];
            [M, merged] = flpm.internal.coalesce(S);
            tc.verifyTrue(merged);
            c = (2*r1 + 2*r2 - 0.1*r2)/2;  h = r1;
            tc.verifyEqual(M.a, sqrt(pi*c*h), 'RelTol', 1e-12);
            tc.verifyEqual(M.s, -r1 + c, 'RelTol', 1e-12);
            tc.verifyEqual(M.ext, 5e-6);
            tc.verifyEqual(M.id, 1);
            % a gap, however small, is interaction, not link-up
            S = state([a1 a2], [0.5e-3 0.5e-3], [0, r1 + r2 + 0.05*r2]);
            [~, merged] = flpm.internal.coalesce(S);
            tc.verifyFalse(merged);
        end

        function linkUpNeedsAlignedPlanes(tc)
            a = 0.05e-3; r = a/sqrt(pi);
            S = state([a a], [0.5e-3 0.5e-3], [0 1.5*r]);
            S.y = [0; 1.9*r];                    % H <= 0.5 (l1 + l2)
            [M, merged] = flpm.internal.coalesce(S);
            tc.verifyTrue(merged);
            tc.verifyEqual(M.y, 0.95*r, 'RelTol', 1e-12);
            S.y = [0; 2.1*r];
            [~, merged] = flpm.internal.coalesce(S);
            tc.verifyFalse(merged);
        end

        function coalescenceIsTransitive(tc)
            a = 0.05e-3; r = a/sqrt(pi); d = 1.9*r;    % 1-2, 2-3 overlap; 1-3 not
            S = state([a a a], [1e-3 1e-3 1e-3], [0 d 2*d]);
            S.dK = [1; 3; 2];
            M = flpm.internal.coalesce(S);
            tc.verifyNumElements(M.a, 1);
            tc.verifyEqual(M.id, 2);
            tc.verifyEqual(M.a, sqrt(pi*(r + d)*r), 'RelTol', 1e-12);
        end

        function pennyCrackStressField(tc)
            ps = @flpm.internal.pennyStress;
            tc.verifyEqual(ps([0 0.5 0.9], [0 0 0]), [-1 -1 -1], 'AbsTol', 1e-12);
            rho = [1.2 2 5];
            tc.verifyEqual(ps(rho, 0*rho), (2/pi)*(1./sqrt(rho.^2 - 1) - asin(1./rho)), 'RelTol', 1e-10);
            z = [0.3 1 3];                         % on the axis (closed form)
            tc.verifyEqual(ps(0*z, z), -(2/pi)*(atan(1./z) - z.*(z.^2 - 1)./(1 + z.^2).^2), 'RelTol', 1e-10);
            % continuity off the plane
            tc.verifyEqual(ps(2, 1e-7), ps(2, 0), 'RelTol', 1e-5);
        end

        function interactionFactorValues(tc)
            r = 1; a = sqrt(pi)*r;
            st = @(s, y) struct('a', [a; a], 's', s(:), 'z', [5; 5], 'y', y(:));
            g = flpm.internal.interactionFactor(st([0 2.1], [0 0]));
            tc.verifyEqual(g, [1 1]' * 1.29, 'RelTol', 0.03);       % gap 0.1 r
            g = flpm.internal.interactionFactor(st([0 2.5], [0 0]));
            tc.verifyEqual(g(1), 1.064, 'RelTol', 0.01);              % gap 0.5 r
            g = flpm.internal.interactionFactor(st([0 0], [0 0.25]));
            tc.verifyEqual(g(1), 0.744, 'RelTol', 0.04);              % coaxial, shielding
            g = flpm.internal.interactionFactor(st([0 30], [0 0]));
            tc.verifyEqual(g, [1; 1]);                                % far apart
            % scale invariance
            S = st([0 2.3], [0 0.4]); S2 = S; S2.a = 1e-4*S.a; S2.s = 1e-4*S.s; S2.z = 1e-4*S.z; S2.y = 1e-4*S.y;
            tc.verifyEqual(flpm.internal.interactionFactor(S2), flpm.internal.interactionFactor(S), 'RelTol', 1e-9);
        end

        function newmanVectorised(tc)
            R = [-1 -0.2 0 0.1 0.5 0.8];
            U = flpm.internal.newmanU(R, 2.73);
            for k = 1:numel(R)
                if R(k) > 0, f = 0.2; else, f = 0.3; end
                tc.verifyEqual(U(k), flpm.internal.newmanU(R(k), 2.73, f), 'RelTol', 1e-12);
            end
        end

        function projectedAreaForms(tc)
            r = 0.1;                                       % sphere
            P = pi * r^2 * [1 1 1];
            tc.verifyEqual(flpm.projectedArea('Projections', P, 'Normal', [1 2 3]), sqrt(pi)*r, 'RelTol', 1e-12);
            e = [0.3 0.1 0.05];                            % aligned ellipsoid
            P = pi * [e(2)*e(3), e(1)*e(3), e(1)*e(2)];
            tc.verifyEqual(flpm.projectedArea('Projections', P), sqrt(P(3)), 'RelTol', 1e-12);
            n = [0.3 -0.5 0.8];
            tc.verifyEqual(flpm.projectedArea('Ellipsoid', e, 'Normal', n), ...
                flpm.projectedArea('Projections', P, 'Normal', n), 'RelTol', 1e-12);
            % rotating the ellipsoid and the normal together changes nothing
            Q = expm([0 -0.4 0.2; 0.4 0 -0.7; -0.2 0.7 0]);
            tc.verifyEqual(flpm.projectedArea('Ellipsoid', e, 'Axes', Q, 'Normal', (Q*n')'), ...
                flpm.projectedArea('Ellipsoid', e, 'Normal', n), 'RelTol', 1e-12);
            % long axis along the load: the smallest silhouette
            tc.verifyEqual(flpm.projectedArea('Ellipsoid', e, 'Normal', [1 0 0]), sqrt(pi*e(2)*e(3)), 'RelTol', 1e-12);
        end

        function residualStressShiftsMeanStress(tc)
            mat = flpm.material('Ti64-AB');
            d = flpm.defectSet(0.06, 0.03, 0);             % near-surface pore
            s0 = flpm.fatigueStrength(d, mat, flpm.options()).stress;
            z0 = flpm.fatigueStrength(d, mat, flpm.options('ResidualStress', @(z) 0*z)).stress;
            tc.verifyEqual(z0, s0, 'AbsTol', 0.06);
            sc = flpm.fatigueStrength(d, mat, flpm.options('ResidualStress', @(z) -200 + 0*z)).stress;
            st = flpm.fatigueStrength(d, mat, flpm.options('ResidualStress', @(z) 200 + 0*z)).stress;
            tc.verifyGreaterThan(sc, s0);                  % compression helps
            tc.verifyLessThan(st, s0);                     % tension hurts
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
S = struct('a', a(:), 'z', z(:), 's', s(:), 'y', zeros(n,1), 'ext', zeros(n,1), ...
    'id', (1:n)', 'dK', zeros(n,1));
end
