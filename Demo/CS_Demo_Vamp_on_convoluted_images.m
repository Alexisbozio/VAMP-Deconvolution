%Deconvolution of a blurred and downsampled image with the Vector
%Approximate Message Passing algorithm with the use of the selective inversion
%method and BM3D denoiser

seed=int64(14753);

addpath(genpath('..'));


profile on

filename='lena.png';
ImIn=double(imread(filename));

%Parameters of the problem
delta=1e-6;                  %variance of noise
height=200;                  %length side of original image
width=height;
downsizing=2;               % downsizing factor
blur_var=1;                  % variance of Gaussian kernel
blur_width=blur_var*6-1;     % width of Gaussian kernel (3.03 deviation from the mean enclose 99% of the pdf)

down_height=height/downsizing;
n=height^2;
m=down_height^2;

%Pick a part of the image
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


%Plot
figure(1);subplot(222);
imy = imagesc(reshape(y, [down_height, down_height]));
colormap gray;axis image;colorbar();pause(0.1);title('Transformed');



%Options
opts.channel_prmts=delta;
opts.x0=x0;
opts.t_max=120;
opts.eps_conv=1e-13;
opts.damp_meas = 0.5;
opts.height=height;
opts.width=width;
r=@channel_awgn;            %Using Selected Inversion
opts.channel=r;


%Run Vamp 
fprintf('o Running V-AMP...\n')
[mse x_hat1] = VAMP_image(y, F, opts);

%%%%%%%%%%%%%%%%%%%%%%%%
%Plot Recovered Signals%
%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
subplot(223);
ima1 = imagesc(reshape(x_hat1, [height, height]));
colormap gray;axis image;colorbar();title('estimate with vamp');

%%%%%%%%%%%%%%%%%%%%%%%%
%Plot mse Trajectories%
%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf;
plot(mse,'.-','Displayname',[denoiser,'-VAMP'])
grid on;
legend(gca,'Show','Location','SouthEast')
xlabel('iteration')
ylabel('mse')

profile viewer; 
profile off
