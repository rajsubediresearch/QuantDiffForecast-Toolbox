function [curvesforecastsens, forecast1] = getensemblemodel(weights1, varargin)
% Mixture-of-trajectories ensemble (variable K):
%   - K = numel(varargin) models, each T x M (same T, same M).
%   - Select exactly M * weights1(i) whole trajectories from model i (sum to M).
%   - Pool and summarize with median and 95% PI at each time.
%
% Usage:
%   [curves, forecast] = getensemblemodel([.5 .3 .2], A1, A2, A3);

    % ---- validate weights and models ----
    K = numel(varargin);
    if K < 1
        error('Provide at least one model matrix.');
    end
    if numel(weights1) ~= K
        error('weights1 length (%d) must equal number of model matrices passed (%d).', numel(weights1), K);
    end
    if any(~isfinite(weights1)) || any(weights1 < 0)
        error('weights1 must be nonnegative and finite.');
    end
    if abs(sum(weights1) - 1) > 1e-8
        error('weights1 must sum to 1 (tolerance 1e-8).');
    end

    mats = varargin;
    [T, M] = size(mats{1});
    for i = 1:K
        if isempty(mats{i}), error('Model %d matrix is empty but weight provided.', i); end
        [Ti, Mi] = size(mats{i});
        if Ti ~= T || Mi ~= M
            error('All model matrices must have the same size T x M.');
        end
    end

    % ---- integer allocation that sums exactly to M ----
    raw = M * weights1(:);
    n_i  = floor(raw);
    rem  = M - sum(n_i);
    if rem > 0
        [~, order] = sort(raw - n_i, 'descend');
        n_i(order(1:rem)) = n_i(order(1:rem)) + 1;  % ensures sum(n_i) == M
    end

    % ---- sample whole trajectories and pool ----
    curvesforecastsens = zeros(T, M);  % preallocate
    colptr = 1;
    for i = 1:K
        if n_i(i) > 0
            if n_i(i) > M, error('Allocation n_i(%d) exceeds M; check weights.', i); end
            idx = randperm(M, n_i(i));                 % sample without replacement
            take = mats{i}(:, idx);
            curvesforecastsens(:, colptr:colptr+n_i(i)-1) = take;
            colptr = colptr + n_i(i);
        end
    end

    % ---- summaries ----
    med = median(curvesforecastsens, 2);
    LB  = quantile(curvesforecastsens', 0.025)';  LB = max(LB, 0);
    UB  = quantile(curvesforecastsens', 0.975)';  UB = max(UB, 0);
    forecast1 = [med, LB, UB];
end
