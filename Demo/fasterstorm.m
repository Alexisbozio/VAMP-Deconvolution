function [ mse, recdensity ] = fasterstorm( y,F,block_width,padheight,x,pa,height )
%FASTERSTORM Summary of this function goes here
%   Detailed explanation goes here

epsilon = 1.5;

   
    L1sum = norm(y,1);
    L2sum = sqrt(L1sum);

   
Fbg = cat(2,F,ones(size(F,1),1));
[m, n] = size(Fbg);
   

c = sum(Fbg,1) ./ max(sum(Fbg,1));
c(end)=0;


cvx_begin
    variable X(n)
    minimize(c*X)
    subject to
	X >= 0
	norm( Fbg * X - y, 2 )  <=(epsilon*L2sum)
cvx_end

z=X(1:end-1);

z = reshape(z, [padheight, padheight]);


zpts   = pkfnd(z,0.01,1);
W=reshape(x, [padheight, padheight]);
[I,J] = ind2sub(size(W),find(W>1e-10));
A=[I,J];
D = pdist2(A,zpts);

[assignment,cost]=munkres(D);

mse=cost .* sqrt(pa) ./ block_width;
recdensity=(block_width^2)*nnz(assignment)/(height^2*pa);


%Plot
figure(42);subplot(223);ima1 = imagesc(reshape(z, [padheight, padheight]));
colormap gray;axis image;pause(0.1);colorbar();title('estimate with strom');

end

