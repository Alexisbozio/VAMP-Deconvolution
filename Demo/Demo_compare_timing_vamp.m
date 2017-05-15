%Deconvolution of a blurred and downsampled toy image (only 0 and 1) with
%the Vector Approximate Message Passing algorithm. Comparison in timing
%between the selected inversion and the use of Woosbury's identity

clear
seed=int64(14753);

addpath(genpath('..'));


%Parameters of the problem
rho=0.03;                    %density
delta=1e-6;                 %variance of noise
height=100;                  %length side of original image
block_width=2;                    %downsizing factor              
blur_var=1;                     %variance of the blur
blur_width=blur_var*6-1;
down_height=height/block_width;
n=height^2;
m=down_height^2;


%Generate Image
x  = generate_image( rho , n );

%Plot
figure(1);subplot(221);imx = imagesc(reshape(x, [height, height]));
colormap gray;axis image;colorbar();pause(0.1);title('Original');

%Generate Super-resolution Matrix F
F = proj_sr( n, block_width, blur_width , blur_var );

% Transformed image y
y = F*x + sqrt(delta) .* randn(m,1) ;

%Plot
figure(1);subplot(222);
imy = imagesc(reshape(y, [down_height, down_height]));
colormap gray;axis image;colorbar();pause(0.1);title('Transformed');


%Options
opts.x0=x;
opts.t_max=30;
opts.eps_conv=1e-13;
opts.delta=delta;
opts.damp_meas = 0.5;
opts.height=height;
f = @prior_01;             
opts.prior_prmts=(rho);
opts.prior=f;
opts.channel_prmts=delta;
r=@channel_awgn;            %Using Selected Inversion
s=@channel_wood;            %Using Woodsbury Identity


%Run Vamp with sel inv
opts.channel=r;
fprintf('o Running V-AMP...\n')
tic;
x_hat1 = vamp(y, F, opts);
toc;
fprintf('o The final MSE is %.4g\n', mean((x_hat1 - x).^2))

%Plot
figure(1);subplot(223);ima1 = imagesc(reshape(x_hat1, [height, height]));
colormap gray;axis image;colorbar();title('estimate with vamp');



%Run Vamp with woodbury identity
opts.channel=s;
fprintf('o Running V-AMPcor...\n')
tic;
x_hat2 = vamp(y, F, opts);
toc;
fprintf('o The final MSE is %.4g\n', mean((x_hat2 - x).^2))

%Plot
figure(1);subplot(224);ima2 = imagesc(reshape(x_hat2, [height, height]));
colormap gray;axis image;colorbar();title('estimate with woodbury');


