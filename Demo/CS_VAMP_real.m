clear
seed=int64(14753);
addpath(genpath('..'));
fwhm = 0.10;
molPhotonRate = 10000;
bgdPhotonRate = 0;
dynamicRange = molPhotonRate - bgdPhotonRate;
pa = 0.166^2 ;

%Parameters of the problem
height=64;                  %length side of original image
block_width=2;                    %downsizing factor              

down_height=height/block_width;
n=height^2;
m=down_height^2;

%Generate Super-resolution Matrix F
F = proj_sr_p( n, block_width,pa ,fwhm);
md=1.4;

 
rho=pa*md/(block_width^2);  

%Generate Image
x  = molPhotonRate*generate_image( rho , n );
molIdx=F*logical(x);

%x(~molIdx)=bgdPhotonRate;

% Transformed image y
y=random('poiss',F*x)/dynamicRange;

x=x/dynamicRange;

delta = diag(max((F*x)/molPhotonRate,1e-6));


%Options
opts.x0=x;
opts.t_max=50;
opts.eps_conv=1e-13;
opts.damp_meas = 0.5;
opts.height=height;
f = @prior_01;  
opts.prior_prmts=rho;
opts.prior=f;
opts.channel_prmts=delta;
r=@channel_awgn;



%Run Vamp with sel inv
opts.channel=r;
fprintf('o Running V-AMP...\n')
tic;
x_hat1 = vamp(y, F, opts);
toc;
fprintf('o The final MSE is %.4g\n', mean((x_hat1 - x).^2))


