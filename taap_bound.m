function lb = taap_bound(m, n, field)
%TAAP_BOUND  Lower bound on mutual coherence: max of six packing bounds.
%  MATLAB port of taap-main/taap_utils.py lowerBound (Massion & Massart).
%  field = 'real' (default) | 'complex'.
%  Bounds: higher-order Welch (1974), Rankin orthoplex (1955), Kabatiansky-
%  Levenshtein (1978), Bukh-Cox (2020), Xia et al. (2005), Bajwa et al. (2012;
%  real only). Stronger than the plain Welch bound (core/getBound) whenever
%  n is large relative to m.
    if nargin < 3 || isempty(field), field = 'real'; end
    isreal_ = strcmpi(field, 'real');
    if m <= 1, error('taap_bound: need m > 1'); end
    if n <= 0, error('taap_bound: need n > 0'); end
    if m >= n, lb = 0; return; end

    % higher-order Welch (degree-k), keep increasing k while it improves
    welch_best = 0;
    deg = 1;
    welch_k = sqrt((n-m) / ((n-1)*m));
    while welch_k > welch_best
        welch_best = welch_k;
        deg = deg + 1;
        binom = exp(gammaln(m+deg) - gammaln(deg+1) - gammaln(m));
        rad = (n/binom - 1) / (n - 1);
        welch_k = max(rad, 0)^(1/(2*deg));
    end

    orthoplex = 0; levenshtein = 0;
    if isreal_ && n > m*(m+1)/2
        orthoplex   = sqrt(1/m);
        levenshtein = sqrt((3*n - m^2 - 2*m) / ((m+2)*(n-m)));
    elseif ~isreal_ && n > m^2
        orthoplex   = sqrt(1/m);
        levenshtein = sqrt((2*n - m^2 - m) / ((m+1)*(n-m)));
    end

    if isreal_
        bukh_cox = (n-m)*(n-m+1) / ...
            (2*n + (n^2 - m*n - n)*sqrt(2+n-m) - (n-m)*(n-m+1));
    else
        bukh_cox = (n-m)^2 / (n + (n^2 - m*n - n)*sqrt(1+n-m) - (n-m)^2);
    end

    xia = 0;
    if log2(n) > m-1, xia = 1 - 2*n^(-1/(m-1)); end

    bajwa = 0;
    if isreal_
        coeff = 2^(2-m)/n / beta(m/2, m/2);
        bajwa = max(0, cos(pi * coeff^(1/(m-1))));
    end

    lb = max([welch_best, orthoplex, levenshtein, bukh_cox, xia, bajwa]);
end
