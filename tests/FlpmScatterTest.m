classdef FlpmScatterTest < matlab.unittest.TestCase
    %FLPMSCATTERTEST  Threshold scatter, copula, truncated fits, bootstrap.

    properties
        mat
        map
    end

    methods (TestClassSetup)
        function buildMap(tc)
            tc.mat = flpm.material('Ti64-AB');
            tc.map = flpm.arrestMap(tc.mat, [], 'SizeRange', [0.01 0.4], 'NumSizes', 16);
        end
    end

    methods (Test)
        function thresholdScatter(tc)
            m = tc.mat; m.dKth_cov = 0.1;
            [mats, w] = flpm.internal.thresholdNodes(m);
            tc.verifyNumElements(mats, 7);
            tc.verifyEqual(sum(w), 1, 'AbsTol', 1e-12);
            f = cellfun(@(q) q.dKth_lc, mats) / m.dKth_lc;
            tc.verifyEqual(sum(w .* log(f).^2), log(1.01), 'RelTol', 1e-9);  % lognormal variance
            % Monte Carlo: per-simulation draw, independent of the sampling stream
            p = flpm.defectPopulation('volumetric', 'Size', 0.06, 'Depth', 0.5, 'Count', 1);
            T0 = flpm.monteCarlo(p, tc.mat, 'SurfaceLength', 10, 'NumSimulations', 40, 'Seed', 2);
            T1 = flpm.monteCarlo(p, m, 'SurfaceLength', 10, 'NumSimulations', 40, 'Seed', 2);
            tc.verifyLessThan(std(T0.stress), 1e-6);
            tc.verifyGreaterThan(std(T1.stress), 5);
            tc.verifyEqual(median(T1.stress), T0.stress(1), 'RelTol', 0.05);
            % closed form: mixture over the 7 nodes
            Ms = cellfun(@(q) flpm.arrestMap(q, [], 'SizeRange', [0.055 0.065], 'NumSizes', 3), ...
                mats, 'UniformOutput', false);
            W = flpm.weakestLink(p, m, 'SurfaceLength', 10, 'Map', Ms);
            tc.verifyEqual(sum(W.killers.weight), 1, 'AbsTol', 1e-9);
            tc.verifyEqual(W.median, T0.stress(1), 'RelTol', 0.03);
            tc.verifyEqual(W.quantile(0.9) - W.quantile(0.1), ...
                quantile(T1.stress, 0.9) - quantile(T1.stress, 0.1), 'RelTol', 0.5);
        end

        function copulaSamplingAndClosedForm(tc)
            base = {'Size', makedist('Lognormal', 'mu', log(0.05), 'sigma', 0.3), ...
                'Depth', makedist('Exponential', 'mu', 0.4), 'Count', 2000, ...
                'SizeBounds', [0.012 0.35]};
            p = flpm.defectPopulation('volumetric', base{:}, 'Copula', {'Gaussian', 0.6});
            s = flpm.samplePopulation(p, 'simulation', 1, 'SurfaceLength', 10, 'Seed', 1);
            tc.verifyEqual(corr(s.size, s.depth, 'Type', 'Spearman'), 6/pi*asin(0.3), 'AbsTol', 0.04);
            pt = flpm.defectPopulation('volumetric', base{:}, 'Copula', {'t', -0.5, 4});
            st = flpm.samplePopulation(pt, 'simulation', 1, 'SurfaceLength', 10, 'Seed', 1);
            tc.verifyLessThan(corr(st.size, st.depth, 'Type', 'Spearman'), -0.35);
            % closed form: rho = 0 equals independence; dependence changes it
            base{6} = 6;
            q0 = flpm.defectPopulation('volumetric', base{:});
            q1 = flpm.defectPopulation('volumetric', base{:}, 'Copula', {'Gaussian', 1e-9});
            q2 = flpm.defectPopulation('volumetric', base{:}, 'Copula', {'Gaussian', -0.8});
            a = {tc.mat, 'SurfaceLength', 10, 'Map', tc.map, 'NumQuad', [120 60 1 1]};
            W0 = flpm.weakestLink(q0, a{:}); W1 = flpm.weakestLink(q1, a{:}); W2 = flpm.weakestLink(q2, a{:});
            tc.verifyEqual(W1.median, W0.median, 'RelTol', 1e-3);
            tc.verifyNotEqual(round(W2.median), round(W0.median));
            tc.verifyError(@() flpm.defectPopulation('surface', 'Size', 0.01, 'Count', 2, ...
                'Copula', {'Gaussian', 0.5}), 'flpm:population:copula');
        end

        function truncatedFit(tc)
            rs = RandStream('mrg32k3a', 'Seed', 4);
            x = exp(log(0.03) + 0.5 * randn(rs, 4000, 1));
            u = 0.03;
            x = x(x > u);
            F = flpm.fitSizeDistribution(x, 'Lognormal', 'DetectionLimit', u);
            tc.verifyEqual(F.dist.mu, log(0.03), 'AbsTol', 0.08);
            tc.verifyEqual(F.dist.sigma, 0.5, 'RelTol', 0.1);
            tc.verifyGreaterThan(F.naive.mu, log(0.03) + 0.2);       % biased
            tc.verifyEqual(F.fractionBelow, 0.5, 'AbsTol', 0.06);
        end

        function dependenceTest(tc)
            rs = RandStream('mrg32k3a', 'Seed', 5);
            z = randn(rs, 300, 2);
            Di = flpm.testDependence(z(:, 1), z(:, 2));
            tc.verifyGreaterThan(Di.spearman(2), 0.01);
            Z = z * chol([1 0.7; 0.7 1]);
            Dd = flpm.testDependence(exp(Z(:, 1)), Z(:, 2));
            tc.verifyTrue(Dd.warranted);
            tc.verifyEqual(Dd.gaussian.rho, 0.7, 'AbsTol', 0.07);
            p = flpm.defectPopulation('volumetric', 'Size', 0.05, 'Depth', 1, 'Count', 2, 'Copula', Dd);
            tc.verifyNotEmpty(p.copula);
        end

        function bootstrapBands(tc)
            rs = RandStream('mrg32k3a', 'Seed', 6);
            data = table(exp(log(0.05) + 0.3 * randn(rs, 80, 1)), 0.05 + 0.6 * rand(rs, 80, 1), ...
                'VariableNames', {'size', 'depth'});
            build = @(t) flpm.defectPopulation('volumetric', 'Size', makedist('Lognormal', ...
                'mu', mean(log(t.size)), 'sigma', std(log(t.size))), 'Depth', ...
                makedist('Uniform', 'lower', 0, 'upper', max(t.depth)), 'Count', 8, ...
                'SizeBounds', [0.012 0.35]);
            B = flpm.bootstrapStrength(data, build, tc.mat, 'SurfaceLength', 10, ...
                'NumBootstrap', 20, 'Seed', 1, 'Map', tc.map);
            tc.verifySize(B.replicates, [20 3]);
            tc.verifyTrue(all(B.ci95(1, :) <= B.estimate + 1 & B.ci95(2, :) >= B.estimate - 1));
            tc.verifyGreaterThan(B.ci95(2, 2) - B.ci95(1, 2), 0);
        end
    end
end
