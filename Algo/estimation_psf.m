%Demo comparing timing between different matrix inversion.
function[mse, eltime, t_final]= estimation_psf(rho,airy,height)

seed=int64(14753);

addpath(genpath('..'));



%Parameters of the problem
delta=1e-6;                 %variance of noise

alpha=2;                    %downsizing factor              

down_height=height/alpha;
n=height^2;
m=down_height^2;


%Generate Image
x  = generate_image( rho , n );


%Generate Super-resolution Matrix F
F = proj_sr( n, alpha, height, airy );

% Transformed image y
y = F*x + sqrt(delta) .* randn(m,1) ;



%Options
opts.x0=x;
opts.t_max=70;
opts.eps_conv=1e-13;
opts.delta=delta;
opts.damp_meas = 0.5;
opts.height=height;
f = @prior_01;             
opts.prior_prmts=(rho);
opts.prior=f;
opts.channel_prmts=delta;
r=@channel_awgn;            %Using Selected Inversion
opts.channel=r;   


%Run Vamp with sel inv


tic;
[x_hat1 t_final] = vamp(y, F, opts);
eltime=toc;
mse= mean((x_hat1 - x).^2);



