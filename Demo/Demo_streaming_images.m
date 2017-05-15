%Super_resolution reconstruction with Vamp on pictures

clear
seed=int64(14753);

addpath(genpath('..'));


profile on

denoiser='BM3D';%Available options are NLM, Gauss, Bilateral, BLS-GSM, BM3D, fast-BM3D, and BM3D-SAPCA 
filename='lena.png';
ImIn=double(imread(filename));

%Parameters of the problem
delta=1e-6;                  %variance of noise
height=40;                  %length side of original image
width=height;
downsizing=2;               % downsizing factor
blur_var=1;                  % variance of Gaussian kernel
blur_width=blur_var*6-1;     % width of Gaussian kernel (3.03 deviation from the mean enclose 99% of the pdf)

down_height=height/downsizing;
n=height^2;
m=down_height^2;

%Create the image
init_1=200;
init_2=200;
x_0=ImIn(init_1:(init_1-1+height),init_2:(init_2-1+width));
Imsize=size(x_0);
x0=x_0(:);
dr = (max(x0)-min(x0));
x0=(x0-min(x0))./dr;

%Rescale the image in range [0,1]
x_0=x_0./dr;      

%Plot
figure(1);subplot(221);imagesc(x_0);title('Original Image');
    axis image;colormap gray;colorbar();pause(0.1);
    
%Generate Super-resolution Matrix F
F= proj_sr( n, downsizing, blur_width, blur_var );

% Transformed image y
y = F*x0 + sqrt(delta) .* randn(m,1) ;

y_square=reshape(y, [down_height, down_height]);

%Plot
figure(1);subplot(222);
imy_sq = imagesc(y_square);
colormap gray;axis image;colorbar();pause(0.1);title('Transformed');


%Block configuration
n_blocks=2;
size_block=down_height/n_blocks;

Bp=sparse(zeros(n,1));
Ap=sparse(zeros(n,1));

for j=1:n_blocks
    for i=1:n_blocks
            
            y_cp=y_square;
            y_cp((i-1)*size_block+1:i*size_block,(j-1)*size_block+1:j*size_block)=0;
            y_block=sparse(y_square-y_cp);
            figure(1);subplot(223);
            imy_block = imagesc(y_block);
            colormap gray;axis image;colorbar();pause(0.1);title('Block');
   
            y_vec=reshape(y_block,[],1);
            
            F_ij=F;
            F_ij(~logical(y_vec),:)=0;
       
            opts.channel_prmts=delta;
            opts.x0=x0;
            opts.t_max=40;
            opts.eps_conv=1e-13;
            opts.damp_meas = 0.5;
            opts.height=height;
            opts.width=width;
            opts.denoiser=denoiser;
            r=@channel_awgn;            %Using Selected Inversion
            s=@channel_wood;            %Using Woodsbury Identity
            opts.channel=r;
            
            
            fprintf('o Running V-AMP...\n')
            tic;
            [x_hat1, A1,B1] = vamp_stream_image(y_vec, F_ij, Ap, Bp, opts);
            Ap=Ap+A1;
            Bp=Bp+B1;
            toc;
            fprintf('o The final MSE is %.4g\n', mean((x_hat1 - x0).^2))
            figure(3);subplot(221);ima = imagesc(reshape(x_hat1, [height, height]));
            colormap gray;axis image;colorbar();title('estimate with vamp');
    end
end    
   
    
          
        


