function [ D ] = generate_down_matrix( height ,alpha )
down_height=height/alpha;
n=height^2;
m=down_height^2;

beta=n/(m);
v = ones(1,alpha)./beta;
K = [];
G=sparse(K);
for i = 1:down_height

    G = blkdiag(G,v);

end

L=sparse(repmat(G,1,alpha));
E=[];
D=sparse(E);
for i = 1:down_height

    D = blkdiag(D,L);

end



end

