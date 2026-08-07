function w = stable_softmax(z)
%STABLE_SOFTMAX  w_i = e^{z_i} / (1 + sum_k e^{z_k}), computed stably.
%  The extra "1 +" is the anchor/self term e^{0}; the shift includes 0 so the
%  result is numerically stable regardless of the magnitude of z. Entries with
%  z_i = -inf (e.g. the excluded anchor atom) map to w_i = 0.
%
%  Ported unchanged from class-idb (src/barriers), which validated it against
%  the per-atom Multi-Similarity / N-pair gradients.
    z  = z(:);
    mx = max([0; z]);
    ez = exp(z - mx);
    ez(~isfinite(z)) = 0;              % z_i = -inf  ->  weight 0
    denom = exp(-mx) + sum(ez);        % = e^{-mx}*(1 + sum e^{z})
    w = ez / denom;
end
