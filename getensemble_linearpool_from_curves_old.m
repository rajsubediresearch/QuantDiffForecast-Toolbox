function [Qens, forecast1, Qmodels, meta] = getensemble_linearpool_from_curves(weights1, varargin)
% Weighted Predictive Density Ensemble (Linear Pool) from per-model curves.
%   - Inputs: K model matrices A1..AK, each T x M (same T, M), plus weights1.
%   - Internally computes per-model quantiles at specified q-levels.
%   - Forms F_ens(x) = sum_i w_i F_i(x) at each time, then inverts to get ensemble quantiles.
%
% Usage:
%   [Qens, forecast1] = getensemble_linearpool_from_curves(w, A1, A2, A3);
%   [Qens, forecast1] = getensemble_linearpool_from_curves(w, A1, A2, 'QLevels', 0.01:0.01:0.99, 'ProbsOut', [0.5 0.025 0.975]);
%
% Outputs:
%   Qens       : T x Lout matrix of ensemble quantiles at 'ProbsOut'
%   forecast1  : T x 3 matrix [median, LB_2.5%, UB_97.5%]
%   Qmodels    : 1xK cell, each T x L matrix of model quantiles used
%   meta       : struct with fields: QLevels, ProbsOut, K, T, M, Tol, MaxIter

    % ---------- split varargin into model matrices vs name-value ----------
    isName = cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), varargin);
    firstNV = find(isName, 1, 'first');
    if isempty(firstNV)
        mats = varargin;
        NV   = {};
    else
        mats = varargin(1:firstNV-1);
        NV   = varargin(firstNV:end);
    end

    % ---------- validate weights & models ----------
    K = numel(mats);
    if K < 1
        error('Provide at least one model matrix (T x M).');
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
    [T, M] = size(mats{1});
    for i = 1:K
        if isempty(mats{i}), error('Model %d matrix is empty but weight provided.', i); end
        [Ti, Mi] = size(mats{i});
        if Ti ~= T || Mi ~= M
            error('All model matrices must have the same size T x M.');
        end
    end
    w = weights1(:)';  % row

    % ---------- parse name-value options ----------
    p = inputParser;
    addParameter(p, 'QLevels', 0.01:0.01:0.99, @(v)isvector(v)&&all(v>0)&all(v<1)&&issorted(v));
    addParameter(p, 'ProbsOut', [0.5 0.025 0.975], @(v)isvector(v)&&all(v>=0)&all(v<=1));
    addParameter(p, 'Tol', 1e-6, @(x)isscalar(x)&&x>0);
    addParameter(p, 'MaxIter', 60, @(x)isscalar(x)&&x>=10);
    addParameter(p, 'ClampNonneg', true, @(x)islogical(x)||ismember(x,[0 1]));
    parse(p, NV{:});
    qlevels = p.Results.QLevels(:)';     % 1 x L
    ProbsOut = p.Results.ProbsOut(:)';   % 1 x Lout
    Tol = p.Results.Tol; MaxIter = p.Results.MaxIter;
    ClampNonneg = p.Results.ClampNonneg;

    L = numel(qlevels);
    Lout = numel(ProbsOut);

    % ---------- derive per-model quantile curves (T x L) ----------
    % Use prctile for broad MATLAB compatibility; quantiles along dim=2 (across columns/samples)
    Qmodels = cell(1, K);
    for i = 1:K
        % Percentile values must be 0..100
        Qi = prctile(mats{i}, 100*qlevels, 2);  % returns T x L
        % Ensure monotone along quantile dimension (guard against numeric jitter)
        Qmodels{i} = enforce_monotone_rows(Qi);
    end

    % ---------- set per-time brackets for inversion ----------
    xmin = inf(T,1); xmax = -inf(T,1);
    for i = 1:K
        xmin = min(xmin, Qmodels{i}(:,1));
        xmax = max(xmax, Qmodels{i}(:,end));
    end
    span = max(xmax - xmin, eps);
    xmin = xmin - 0.05*span;  % small padding
    xmax = xmax + 0.05*span;

    % ---------- compute ensemble quantiles at ProbsOut ----------
    Qens = nan(T, Lout);
    for t = 1:T
        % cache knots for all models at time t
        xk = cell(1,K);
        pk = cell(1,K);
        for i = 1:K
            xk{i} = Qmodels{i}(t, :);   % 1 x L, nondecreasing
            pk{i} = qlevels;            % 1 x L, strictly increasing
        end
        for a = 1:Lout
            p_out = ProbsOut(a);
            Qens(t,a) = invert_pooled_cdf(p_out, xk, pk, w, xmin(t), xmax(t), Tol, MaxIter);
        end
    end

    % ---------- default summary [median, LB, UB] ----------
    % find indices if ProbsOut already includes desired probs
    p_needed = [0.5, 0.025, 0.975];
    pos = zeros(1,3);
    for j = 1:3
        ix = find(abs(ProbsOut - p_needed(j)) < 1e-12, 1);
        pos(j) = iff(isempty(ix), 0, ix);
    end
    forecast1 = nan(T,3);
    for t = 1:T
        for j = 1:3
            if pos(j) > 0
                forecast1(t,j) = Qens(t, pos(j));
            else
                % compute on the fly at this single prob
                forecast1(t,j) = invert_pooled_cdf(p_needed(j), ...
                    arrayfun(@(i) Qmodels{i}(t,:), 1:K, 'uni', false), ...
                    repmat({qlevels}, 1, K), w, xmin(t), xmax(t), Tol, MaxIter);
            end
        end
    end
    if ClampNonneg
        forecast1(:,2) = max(forecast1(:,2), 0);
        forecast1(:,3) = max(forecast1(:,3), 0);
    end

    if nargout >= 4
        meta = struct('QLevels', qlevels, 'ProbsOut', ProbsOut, ...
                      'K', K, 'T', T, 'M', M, 'Tol', Tol, 'MaxIter', MaxIter);
    end
end

% ===================== helpers =====================

function x = invert_pooled_cdf(p_out, xk_cell, pk_cell, w, xl, xr, tol, maxiter)
% Find x in [xl, xr] such that F_mix(x) = p_out, where F_mix = sum_i w_i F_i(x).
% Each F_i is piecewise-linear CDF derived from per-model (x_k, p_k) knots.

    fl = pooled_cdf(xl, xk_cell, pk_cell, w) - p_out;
    fr = pooled_cdf(xr, xk_cell, pk_cell, w) - p_out;
    if fl > 0 && fr > 0, x = xl; return; end   % below left tail
    if fl < 0 && fr < 0, x = xr; return; end   % above right tail

    % Bisection
    for it = 1:maxiter
        xm = 0.5*(xl + xr);
        fm = pooled_cdf(xm, xk_cell, pk_cell, w) - p_out;
        if abs(fm) < tol || (xr - xl) < tol*max(1,abs(xm))
            x = xm; return;
        end
        if sign(fm) == sign(fl)
            xl = xm; fl = fm;
        else
            xr = xm; fr = fm;
        end
    end
    x = 0.5*(xl + xr); % fallback
end

function F = pooled_cdf(x, xk_cell, pk_cell, w)
% F(x) = sum_i w_i * F_i(x), with F_i(x) piecewise-linear between (x_k, p_k)
    K = numel(w);
    Fi = zeros(1,K);
    for i = 1:K
        Fi(i) = model_cdf_from_quantiles(x, xk_cell{i}, pk_cell{i});
    end
    F = sum(w .* Fi);
end

function p = model_cdf_from_quantiles(x, xk, pk)
% Given knots (xk, pk) with xk nondecreasing and pk strictly increasing,
% return F(x) by linear interpolation in (x,p). Flat tails outside knots.
    L = numel(xk);
    if x <= xk(1)
        p = pk(1);
        return;
    elseif x >= xk(L)
        p = pk(L);
        return;
    else
        j = find(xk <= x, 1, 'last');
        if j == L, p = pk(L); return; end
        x0 = xk(j); x1 = xk(j+1); p0 = pk(j); p1 = pk(j+1);
        if x1 == x0
            p = p1; return;
        end
        t = (x - x0) / (x1 - x0);
        p = (1 - t)*p0 + t*p1;
    end
end

function Qmono = enforce_monotone_rows(Q)
% Ensure nondecreasing quantiles along columns for each row using cumulative max.
    Qmono = Q;
    for r = 1:size(Q,1)
        Qmono(r,:) = cummax(Q(r,:));
    end
end

function y = iff(cond, a, b)
% inline conditional
    if cond, y = a; else, y = b; end
end
