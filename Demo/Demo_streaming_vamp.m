%Deconvolution of a blurred and downsampled toy image (only 0 and 1) with
%the Vector Approximate Message Passing algorithm using the woodbury
%identity method and by dividing the image by block

clear
seed=int64(14753);

addpath(genpath('..'));


%Parameters of the problem
rho=0.03;                    %density
delta=1e-6;                 %variance of noise
height=100;                  %length side of original image
downsizing=2;                    %downsizing factor              
blur_var=1;                     %variance of the blur
blur_width=blur_var*6-1;
down_height=height/downsizing;
n=height^2;
m=down_height^2;


%Generate Image
x  = generate_image( rho , n );

%Plot
figure(1);subplot(221);imx = imagesc(reshape(x, [height, height]));
colormap gray;axis image;colorbar();pause(0.1);title('Original');

%Generate Super-resolution Matrix F
F = proj_sr(n, downsizing, blur_width , blur_var);

% Transformed image y
y = F*x + sqrt(delta) .* randn(m,1) ;
y_square=reshape(y, [down_height, down_height]);

%Plot
figure(1);subplot(222);
imy_sq = imagesc(y_square);
colormap gray;axis image;colorbar();pause(0.1);title('Transformed');


%Block configuration
n_blocks=5;
size_block=down_height/n_blocks;



%Vamp options

opts.t_max=30;
opts.eps_conv=1e-13;
opts.delta=delta;
opts.damp_meas = 0.5;
f = @prior_01;  
opts.prior_prmts=rho;
opts.prior=f;
opts.channel_prmts=delta;
r=@channel_awgn;
%opts.channel=r;
s=@channel_wood;            %Using Woodsbury Identity
opts.channel=s;

x_hat1=zeros(n,1);
A=zeros(n,1);
B=zeros(n,1);

for j=1:n_blocks
    for i=1:n_blocks
            
            
            
            %Select the block
            rows=((i-1)*size_block+1:i*size_block);
            column=((j-1)*size_block+1:j*size_block);
            block=zeros(down_height);
            block(rows,column)=1;
            sel_el_m=logical(reshape(block,[],1));
            
          
            y_vec=y(sel_el_m,:);
            
            
            %Build the matrix that is talking to the block
            F_ij=F(sel_el_m,:);
            sel_el_n=logical(sum(logical(F_ij)));
            
            m_red=nnz(sel_el_m);
            n_red=nnz(sel_el_n);
                    
            F_red = F_ij(:, sel_el_n);
            Ap=A(sel_el_n,:);
            Bp=B(sel_el_n,:);
            
            x0=x(sel_el_n);
            opts.x0=x0;
            %Running Vamp
            fprintf('o Running V-AMP...\n')
            tic;
            [x_hat_red, A1, B1] = vamp_stream_wood(y_vec, F_red, Ap, Bp, opts);
            toc;
            
            
            %Update priors
            A(sel_el_n)=A(sel_el_n)+A1;
            B(sel_el_n)=B(sel_el_n)+B1;
            x_hat1(sel_el_n)=x_hat_red;
            
                              
            %Print
            figure(3);subplot(221);ima = imagesc(reshape(x_hat1, [height, height]));
            colormap gray;axis image;colorbar();title('estimate with vamp');
    end
end    

    
fprintf('o The final MSE is %.4g\n', mean((x_hat1 - x).^2))
        


