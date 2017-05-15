%Super_resolution reconstruction with Vamp and the
%Generalized Expectation Consistent Signal Recovery for Nonlinear
%Measurements.

clear
seed=int64(14753);

addpath(genpath('..'));


%Parameters of the problem
rho=0.04;                   %density
delta=1e-6;                 %variance of noise
height=30;                  %length side of original image
alpha=2;                    %downsizing factor              
blur_var=1;                     %size of the PSF
blur_width=blur_var*6-1; 
down_height=height/alpha;
n=height^2;
m=down_height^2;


%Generate Image
x  = generate_image( rho , n );

%Plot
figure(1);subplot(221);imx = imagesc(reshape(x, [height, height]));
colormap gray;axis image;colorbar();pause(0.1);title('Original');

%Generate Super-resolution Matrix F
C = generate_convolution_matrix(  n, blur_width, blur_var );        %Blurring operator 
D = generate_down_matrix( height , alpha );       %Downsizing Matrix
F = D*C;

% Transformed image y
y = F*x + sqrt(delta) .* randn(m,1) ;
%Plot
figure(1);subplot(222);
imy = imagesc(reshape(y, [down_height, down_height]));
colormap gray;axis image;colorbar();pause(0.1);title('Transformed');


%Options
opts.x0=x;
opts.t_max=50;
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

opts.channel=r;
%opts.channel=s;


%Run Vamp 
fprintf('o Running V-AMP...\n')

x_hat1 = vamp(y, F, opts);

fprintf('o The final MSE is %.4g\n', mean((x_hat1 - x).^2))

%Plot
figure(1);subplot(223);ima1 = imagesc(reshape(x_hat1, [height, height]));
colormap gray;axis image;colorbar();title('estimate with vamp');


%Run Multi Layer Vamp of the Corean type
fprintf('o Running V-AMPcor...\n')

x_hat2 = vamp_corean(y, C, D, opts);

fprintf('o The final MSE is %.4g\n', mean((x_hat2 - x).^2))

%Plot
figure(1);subplot(224);ima2 = imagesc(reshape(x_hat2, [height, height]));
colormap gray;axis image;colorbar();title('estimate with multi layer vamp');


