function [a,c]= prior_gb(A, B, prior_prmts)
    [rho, mu, sig]=deal(prior_prmts{:});
    m = (B * sig + mu) ./ (1. + A * sig);
    v = sig ./ (1 + A * sig);
    p_s = rho ./ ( rho + (1 - rho) * sqrt(1. + A .* sig) .*exp(-.5 * m.^2 ./ v + .5 * mu.^2 ./ sig) );

    a = p_s .* m;
    c = p_s .* v + p_s .* (1. - p_s) .* m.^2;
    c = max(c,1e-10);
  
end