classdef FlpmValidationTest < matlab.unittest.TestCase
    %FLPMVALIDATIONTEST  S-N maximum likelihood and killer-defect statistics.

    methods (Test)
        function rflRecoversParameters(tc)
            [S, N, ro] = simulateRFL(7, 400);
            F = flpm.fitSN(S, N, ro);
            tc.verifyEqual(F.params.b1, -2.5, 'AbsTol', 0.5);
            tc.verifyEqual(exp(F.params.mu_g), 250, 'RelTol', 0.05);
            tc.verifyEqual(F.params.s_g, 0.08, 'AbsTol', 0.04);
            tc.verifyEqual(F.params.s, 0.3, 'RelTol', 0.25);
            tc.verifySize(F.se, [1 5]);
            % fitted life quantiles are consistent with the cdf
            n50 = F.quantile(320, 0.5);
            tc.verifyEqual(F.cdf(320, n50), 0.5, 'AbsTol', 1e-6);
        end

        function modelLimitComparison(tc)
            [S, N, ro] = simulateRFL(8, 300);
            rs = RandStream('mrg32k3a', 'Seed', 1);
            right = 250 * exp(0.08 * randn(rs, 5000, 1));      % true limit
            wrong = 200 * exp(0.08 * randn(rs, 5000, 1));      % biased
            Ff = flpm.fitSN(S, N, ro);
            C1 = flpm.compareSN(Ff, flpm.fitSN(S, N, ro, 'FatigueLimit', right));
            C2 = flpm.compareSN(Ff, flpm.fitSN(S, N, ro, 'FatigueLimit', wrong));
            tc.verifyGreaterThan(C1.p, 0.01);
            tc.verifyLessThan(C2.p, 1e-3);
            tc.verifyGreaterThan(C1.deltaAIC, C2.deltaAIC);
        end

        function killerKsTest(tc)
            T = table((1:400)', 0.05 + 0.02 * rand(RandStream('mrg32k3a', 'Seed', 2), 400, 1), ...
                linspace(0, 1, 400)', 'VariableNames', {'simulation', 'killerSize', 'killerDepth'});
            rs = RandStream('mrg32k3a', 'Seed', 3);
            same = [0.05 + 0.02 * rand(rs, 30, 1), rand(rs, 30, 1)];
            off = [0.08 + 0.02 * rand(rs, 30, 1), 0.5 * rand(rs, 30, 1)];
            K1 = flpm.killerDistribution(T, 'Observed', same);
            K2 = flpm.killerDistribution(T, 'Observed', off);
            tc.verifyGreaterThan(K1.ksSize(2), 0.05);
            tc.verifyGreaterThan(K1.ksDepth(2), 0.05);
            tc.verifyLessThan(K2.ksSize(2), 1e-6);
            tc.verifyLessThan(K2.ksDepth(2), 0.01);
            tc.verifyEqual(K1.sizeQuantiles(3), 0.06, 'AbsTol', 0.003);
        end
    end
end

function [S, N, ro] = simulateRFL(seed, n)
rs = RandStream('mrg32k3a', 'Seed', seed);
levels = [260 280 300 330 370 420];
S = levels(randi(rs, numel(levels), n, 1))';
g = 250 * exp(0.08 * randn(rs, n, 1));
N = exp(22 - 2.5 * log(max(S - g, eps)) + 0.3 * randn(rs, n, 1));
N(S <= g) = Inf;
ro = N > 1e7;
N(ro) = 1e7;
end
