function [ x ] = generate_image( rho , n )

k=int64(ceil(rho * n));     %k=sparsity
x = zeros(n,1);             %size of the image
beq = randperm(n);          %distributes randomly positions of non zero component in x
supp = beq(1:k);            %build a vector with k first elements of beq (bequille)
x(supp) = 1.;               %put 1 in every position given by supp 

end

