%Deconvolution of a blurred and downsampled toy image (only 0 and 1) with
%the Vector Approximate Message Passing algorithm. Comparison in timing
%between the selected inversion and the use of Woosbury's identity

clear
seed=int64(14753);

addpath(genpath('..'));
fwhm = 0.10;
molPhotonRate = 750;
bgdPhotonRate = 25;
dynamicRange = molPhotonRate - bgdPhotonRate;
pa = 0.166^2 ;
md = 3      ;

%Parameters of the problem
height=60;                      %length side of original image
block_width=2;                  %downsizing factor              
down_height=height/block_width;
n=height^2;
m=down_height^2;

rho=pa*md/(block_width^2);  

vampGBMean = 1.0;
vampGBVar = molPhotonRate ./ (dynamicRange.^2);
vampGBRho = rho;

%Generate Image
x  = molPhotonRate.*generate_image( rho , n );
molIdx=logical(x);
x(~molIdx)=bgdPhotonRate;


%Generate Super-resolution Matrix F
F = proj_sr_p( n, block_width,pa ,fwhm);

% Transformed image y
highFilterVal = max(F*ones(n,1));
z=random('poiss',F*x);
y=(z - bgdPhotonRate*highFilterVal) ./ dynamicRange;

%delta = (F*x) ./ (dynamicRange.^2);
delta = bgdPhotonRate ./ (dynamicRange.^2);    
delta = max(delta,1e-6);

x=x./dynamicRange;

%Plot
figure(1);subplot(221);imx = imagesc(reshape(x, [height, height]));
colormap gray;axis image;colorbar();pause(0.1);title('Original');
%Plot
figure(1);subplot(222);
imy = imagesc(reshape(y, [down_height, down_height]));
colormap gray;axis image;colorbar();pause(0.1);title('Transformed');


%Options
opts.x0=x;
opts.t_max=30;
opts.eps_conv=1e-13;
opts.damp_meas = 0.5;
opts.height=height;
f = @prior_01; 
g= @prior_gb;
opts.prior_prmts=[vampGBRho,vampGBMean,vampGBVar];
opts.prior=g;
opts.channel_prmts=diag(delta);
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



%{
Run Vamp with woodbury identity
opts.channel=s;
fprintf('o Running V-AMPcor...\n')
tic;
x_hat2 = vamp(y, F, opts);
toc;
fprintf('o The final MSE is %.4g\n', mean((x_hat2 - x).^2))

%Plot
figure(1);subplot(224);ima2 = imagesc(reshape(x_hat2, [height, height]));
colormap gray;axis image;colorbar();title('estimate with woodbury');

%}
