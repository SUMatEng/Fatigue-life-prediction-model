classdef FlpmPhysicsTest < matlab.unittest.TestCase
    %FLPMPHYSICSTEST  Unit tests of the constitutive pieces of the model.
    %   Run with:  results = runtests('tests');

    methods (Test)
        function newmanClosureRatio(tc)
            % Value of the published closure function for Ti64-AB, R = 0.1.
            U = flpm.internal.newmanU(0.1, 2.73, []);
            tc.verifyEqual(U, 0.782309004496420, 'AbsTol', 1e-12);
            % Fully open for R close to 1 is not required, but U <= 1 always.
            tc.verifyLessThanOrEqual(flpm.internal.newmanU(0.7, 2.73, []), 1);
        end

        function formanSmallCrackLimit(tc)
            % a -> 0 recovers the semi-circular surface crack value ~0.657.
            Y = flpm.internal.formanShaftY([1e-9 1e-6], 7.3e-3);
            tc.verifyEqual(Y, 0.657144393828712*[1 1], 'RelTol', 1e-3);
            % Monotonic increase and hold beyond D/2.
            a = linspace(1e-5, 4e-3, 50);
            Ya = flpm.internal.formanShaftY(a, 7.3e-3);
            tc.verifyTrue(all(diff(Ya(a < 3.65e-3)) > 0));
            tc.verifyEqual(Ya(end), flpm.internal.formanShaftY(3.65e-3, 7.3e-3), 'RelTol', 1e-12);
        end

        function netSectionStressLimits(tc)
            D = 7.3e-3;
            tc.verifyEqual(flpm.internal.netSectionStress(1e-7, 100, D, 0.13), 100, 'RelTol', 1e-3);
            a = linspace(1e-5, 3e-3, 40);
            S = arrayfun(@(x) flpm.internal.netSectionStress(x, 100, D, 0.13), a);
            tc.verifyTrue(all(diff(S) > 0), 'net-section stress must increase with crack depth');
            % Linear in the remote stress.
            tc.verifyEqual(flpm.internal.netSectionStress(1e-3, 300, D, 0.13), ...
                3*flpm.internal.netSectionStress(1e-3, 100, D, 0.13), 'RelTol', 1e-12);
        end

        function roughnessShieldingFactor(tc)
            [YR, sa] = flpm.roughnessShielding(0.010, 0.100);   % a/2b = 0.1
            tc.verifyEqual(YR, 1.032170528546519, 'AbsTol', 1e-12);
            tc.verifyEqual(sa, 0.1*(2.97*0.1 - 3.51*0.01 - 9.74*0.001), 'AbsTol', 1e-15);
            % Steep valleys: sqrt(areaR) = 0.38*2b.
            [~, sa2] = flpm.roughnessShielding(0.05, 0.1);
            tc.verifyEqual(sa2, 0.038, 'AbsTol', 1e-15);
            % Vector Rv with scalar RSm.
            tc.verifySize(flpm.roughnessShielding([0.01 0.02 0.05], 0.1), [1 3]);
        end

        function classificationAndInternalDamping(tc)
            D = 7.3e-3; o = flpm.modelOptions();
            S = state([0.05 0.05 0.05], [0 3.58 3.62], [0 pi/2 pi], [1 1 1]);
            [Y, cls] = flpm.internal.geometryFactors(S, D, o, 1);
            tc.verifyEqual(cls(:).', [0 1 2]);        % internal, near, surface
            % Surface: Forman bounded below by the surface floor 1.0335.
            tc.verifyEqual(Y(3), 1.0335, 'AbsTol', 1e-12);
            % Internal Y carries the 0.5/0.65 factor relative to the raw value.
            Yraw = flpm.internal.alegreEmbeddedY(S.a(1), 1, D/2 + S.a(1), D/2);
            tc.verifyEqual(Y(1), 0.5/0.65*Yraw, 'RelTol', 1e-12);
        end

        function coalescenceMergesByEquivalentArea(tc)
            % Two defects in contact (centre distance < a1 + a2).
            S = state([0.06 0.04], [1.0 1.09], [0 0], [0.5 1]);
            S.ext = [2e-6; 5e-6];
            S.dK = [2; 1];
            [M, merged] = flpm.internal.coalesce(S, 1);
            tc.verifyTrue(merged);
            tc.verifyNumElements(M.a, 1);
            tc.verifyEqual(M.a, hypot(0.06, 0.04)/1000, 'RelTol', 1e-12);
            tc.verifyEqual(M.x, (0.06*1.0 + 0.04*1.09)/0.10/1000, 'RelTol', 1e-12);
            tc.verifyEqual(M.asp, (0.06*0.5 + 0.04*1)/0.10, 'RelTol', 1e-12);
            tc.verifyEqual(M.ext, 5e-6);
            tc.verifyEqual(M.id, 1);                     % member with max dK
        end

        function coalescenceIsTransitive(tc)
            % A touches B, B touches C, A does not touch C -> one defect.
            S = state([0.05 0.05 0.05], [0 0.09 0.18], [0 0 0], [1 1 1]);
            S.dK = [1; 3; 2];
            M = flpm.internal.coalesce(S, 1);
            tc.verifyNumElements(M.a, 1);
            tc.verifyEqual(M.id, 2);
            tc.verifyEqual(M.a, sqrt(3)*0.05e-3, 'RelTol', 1e-12);
        end

        function noCoalescenceWhenApart(tc)
            S = state([0.05 0.05], [0 0.2], [0 0], [1 1]);
            [M, merged] = flpm.internal.coalesce(S, 1);
            tc.verifyFalse(merged);
            tc.verifyNumElements(M.a, 2);
        end
    end
end

function S = state(sizeMm, xMm, ang, asp)
% Defects placed along the x axis at the given radial distances.
geom = flpm.geometry();
d = flpm.defectSet(sizeMm, abs(xMm), ang, asp);
S = flpm.internal.prepareState(d, geom);
S.dK = zeros(size(S.a));
end
