function [ mse,recdensity ] = vampoisson( y,F,block_width,padheight,x,pa,height,dynamicRange,n,down_height,rho)

x=x/dynamicRange;

delta = max((F*x)/dynamicRange,1/dynamicRange^2);


y=y/dynamicRange;

%Block configuration
n_blocks=1;
size_block=down_height/n_blocks;



%Vamp options
%opts.x0=x;
opts.t_max=70;
opts.eps_conv=1e-13;
opts.damp_meas = 0.5;
f = @prior_01;  
opts.prior_prmts=rho;
opts.prior=f;
%opts.channel_prmts=diag(delta);
r=@channel_awgn;
opts.channel=r;
%s=@channel_wood; 
opts.height=padheight;
%Using Woodsbury Identity
%opts.channel=s;



%x_hat1 = vamp(y, F, opts);





x_hat1=zeros(n,1);
A=zeros(n,1);
B=zeros(n,1);





j_axis=randperm(n_blocks);
i_axis=randperm(n_blocks);

  tic;
  

    for j=j_axis
        for i=i_axis



                %Select the block
                rows=((i-1)*size_block+1:i*size_block);
                column=((j-1)*size_block+1:j*size_block);
                block=zeros(down_height);
                block(rows,column)=1;
                sel_el_m=logical(reshape(block,[],1));
                y_vec=y(sel_el_m,:);


                %Build the reduced matrix F_red that is talking to the block

                F_ij=F(sel_el_m,:);
                opts.channel_prmts=diag(delta(sel_el_m));
                sel_el_n=logical(sum(logical(F_ij)));
                m_red=nnz(sel_el_m);
                n_red=nnz(sel_el_n);
                F_red = F_ij(:, sel_el_n);
                Ap=A(sel_el_n,:);
                Bp=B(sel_el_n,:);
              

                x0=x(sel_el_n,:);
                opts.x0=x0;
                %Running Vamp
                fprintf('o Running V-AMP...\n')
              
                [x_hat_red, A1, B1] = vamp_stream_wood(y_vec, F_red, Ap, Bp, opts);
            


                %Update priors
                A(sel_el_n)=Ap+A1;
                B(sel_el_n)=Bp+B1;
                x_hat1(sel_el_n)=x_hat_red;


                
        end
    end    

    toc;
   
        if mean((x_hat1 - x).^2)==0

            mse=0;
            recdensity=(block_width^2)*nnz(x_hat1)/(height^2*pa);
  

        else    
            z = reshape(x_hat1, [padheight, padheight]);


            zpts   = pkfnd(z,0.01,1);
            W=reshape(x, [padheight, padheight]);
            [I,J] = ind2sub(size(W),find(W>1e-10));
            A=[I,J];
            D = pdist2(A,zpts);

            [assignment,cost]=munkres(D);

            mse=cost .* sqrt(pa) ./ block_width;
            recdensity=(block_width^2)*nnz(assignment)/(height^2*pa);

        end
figure(42);subplot(224);ima1 = imagesc(reshape(x_hat1, [padheight, padheight]));
colormap gray;axis image;pause(0.1);colorbar();title('estimate with vamp');
end

