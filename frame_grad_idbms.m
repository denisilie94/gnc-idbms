function g = frame_grad_idbms(v, d, F, j, p)
%FRAME_GRAD_IDBMS  Atom gradient of the GNC-IDBMS objective.
%  Keeps IDB's weighted quadratic but replaces the hard hinge by the smooth
%  multi-similarity (log-sum-exp) barrier
%
%    B_beta(d) = (1/beta) * log( 1 + sum_{i~=j} exp( beta*(|v_i| - mu) ) ),
%
%  whose gradient is a softmax-weighted repulsion (the 1/beta factor cancels):
%
%    g = F * (max(|v|/mu, 1) .* v)                           (weighted quadratic)
%        + lambda * F * ( omega .* sign(v) ),                (smooth barrier)
%    omega_i = e^{beta(|v_i| - mu)} / ( 1 + sum_k e^{beta(|v_k| - mu)} ).
%
%  Unlike the hinge, every omega_i is strictly positive, so the barrier keeps
%  separating a pair even after it satisfies the target; and the weights are
%  self-paced, concentrating on the most correlated pairs. The temperature beta
%  is annealed across the sweeps by the caller (see idb_frame).
%
%  Input:
%    v  - correlations F'*d of atom j with the others, with v(j) = 0  (n x 1)
%    d  - the current atom F(:,j), unit norm                          (m x 1)
%    F  - the full frame                                              (m x n)
%    j  - index of the current atom
%    p  - parameters: p.mu (target), p.lambda (barrier weight),
%         p.beta (temperature)
%  Output:
%    g  - gradient for atom j                                         (m x 1)

    s = abs(v);

    % weighted quadratic over all pairs (v(j) = 0, so atom j contributes 0)
    g = F*(max(s/p.mu, 1) .* v);

    % smooth multi-similarity barrier
    z = p.beta*(s - p.mu);
    z(j) = -inf;                       % exclude the atom itself
    wb = stable_softmax(z);
    g = g + p.lambda * (F*(wb .* sign(v)));
end
