function g = frame_grad_idb(v, d, F, j, p)
%FRAME_GRAD_IDB  Atom gradient of the distance-barrier objective (IDB).
%  Weighted quadratic plus hinge barrier:
%
%    g = F * (w .* v)                                        (weighted quadratic)
%        + lambda * [ sum_{v_i >  mu} (d_i - d)
%                   - sum_{v_i < -mu} (d_i + d) ]            (hinge barrier)
%    w_i = max(|v_i|/mu, 1)
%
%  The hinge term is active only on the pairs that violate the target, i.e.
%  those with |v_i| > mu; the quadratic term acts on all pairs, since w_i >= 1.
%
%  Input:
%    v  - correlations F'*d of atom j with the others, with v(j) = 0  (n x 1)
%    d  - the current atom F(:,j), unit norm                          (m x 1)
%    F  - the full frame                                              (m x n)
%    j  - index of the current atom
%    p  - parameters: p.mu (target coherence), p.lambda (barrier weight)
%  Output:
%    g  - gradient for atom j                                         (m x 1)

    w = max(abs(v)/p.mu, 1);
    g = F*(w.*v);

    im = v >  p.mu;  nm = nnz(im);      % pairs too close, positive correlation
    ip = v < -p.mu;  np = nnz(ip);      % pairs too close, negative correlation
    if nm
        g = g + p.lambda*(sum(F(:,im),2) - nm*d);
    end
    if np
        g = g - p.lambda*(sum(F(:,ip),2) + np*d);
    end
end
