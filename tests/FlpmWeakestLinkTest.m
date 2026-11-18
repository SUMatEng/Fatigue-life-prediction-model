classdef FlpmWeakestLinkTest < matlab.unittest.TestCase
    %FLPMWEAKESTLINKTEST  Closed-form weakest-link strength distribution.

    properties
        mat
        map
    end

    methods (TestClassSetup)
        function buildMap(tc)
            tc.mat = flpm.material('Ti64-AB');
            tc.map = flpm.arrestMap(tc.mat, [], 'SizeRange', [0.015 0.25], 'NumSizes', 14);
        end
    end

    methods (Test)
        function mapMatchesSolver(tc)
            a = [0.031 0.047 0.083 0.12];
            z = [0.6 0.031 0.06 0.5];         % internal, surface, near-surface, internal
            for k = 1:numel(a)
                r = flpm.fatigueStrength(flpm.defectSet(a(k), z(k), 0), tc.mat);
                tc.verifyEqual(tc.map.eval(a(k), z(k)), r.stress, 'RelTol', 0.02);
            end
            tc.verifyError(@() tc.map.eval(1, 1), 'flpm:arrestMap:range');
        end

        function constantDefectGivesStep(tc)
            p = flpm.defectPopulation('volumetric', 'Size', 0.06, 'Depth', 0.5, 'Count', 3);
            W = flpm.weakestLink(p, tc.mat, 'SurfaceLength', 10, 'Map', tc.map);
            sw = tc.map.eval(0.06, 0.5);
            tc.verifyEqual(W.median, sw, 'RelTol', 1e-12);
            tc.verifyEqual(W.PfFun([0.999 1.001] * sw), [0 1]);
            tc.verifyEqual(sum(W.killers.weight), 1, 'AbsTol', 1e-12);
        end

        function poissonClosedForm(tc)
            % one size class, Poisson count: Pf = 1 - exp(-mu) above sigma_w
            p = flpm.defectPopulation('volumetric', 'Size', 0.06, 'Depth', 0.5, ...
                'Count', 2, 'CountModel', 'poisson');
            W = flpm.weakestLink(p, tc.mat, 'SurfaceLength', 10, 'Map', tc.map);
            tc.verifyEqual(W.Pf(end), 1 - exp(-2), 'RelTol', 1e-12);
            tc.verifyEqual(sum(W.killers.weight), 1 - exp(-2), 'RelTol', 1e-12);
        end

        function stressFieldScales(tc)
            p = flpm.defectPopulation('volumetric', ...
                'Size', makedist('Lognormal', 'mu', log(0.05), 'sigma', 0.3), ...
                'Depth', makedist('Exponential', 'mu', 0.4), 'Count', 6, ...
                'SizeBounds', [0.016 0.24]);
            W1 = flpm.weakestLink(p, tc.mat, 'SurfaceLength', 10, 'Map', tc.map);
            W2 = flpm.weakestLink(p, tc.mat, 'SurfaceLength', 10, 'Map', tc.map, ...
                'StressField', @(z, s, y) 1.25 + 0 * z, 'NumQuad', [200 80 2 1]);
            tc.verifyEqual(W2.median, W1.median / 1.25, 'RelTol', 1e-9);
            % a gradient that unloads deep defects raises the strength
            W3 = flpm.weakestLink(p, tc.mat, 'SurfaceLength', 10, 'Map', tc.map, ...
                'StressField', @(z, s, y) exp(-z), 'NumQuad', [200 80 2 1]);
            tc.verifyGreaterThan(W3.median, W1.median);
        end

        function matchesMonteCarloWithoutInteraction(tc)
            % defects spread so thinly that none interact
            p = flpm.defectPopulation('volumetric', ...
                'Size', makedist('Lognormal', 'mu', log(0.05), 'sigma', 0.3), ...
                'Depth', makedist('Exponential', 'mu', 0.4), 'Count', 5, ...
                'CountModel', 'poisson', 'SizeBounds', [0.016 0.24]);
            W = flpm.weakestLink(p, tc.mat, 'SurfaceLength', 1e4, 'AxialLength', 1e4, ...
                'Map', tc.map);
            T = flpm.monteCarlo(p, tc.mat, 'SurfaceLength', 1e4, 'AxialLength', 1e4, ...
                'NumSimulations', 300, 'Seed', 3, 'UseParallel', true);
            q = [0.1 0.25 0.5 0.75 0.9];
            tc.verifyEqual(W.quantile(q), quantile(T.stress, q), 'RelTol', 0.03);
        end
    end
end
