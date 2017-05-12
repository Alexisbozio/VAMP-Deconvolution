
  % prior_01(A, B, prior_prmts)

 %  The binary prior. 

 function[a,v]= prior_01(A, B, prior_prmts)
    rho = prior_prmts(1);
    a = 1. ./ (1. + exp(-( log(rho ./ (1. - rho)) + B - .5 * A )));
    v=max(1e-12, a - a.^2);
end
