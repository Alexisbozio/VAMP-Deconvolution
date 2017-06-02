clear
seed=int64(14753);
addpath(genpath('..'));
fwhm = 0.10;
molPhotonRate = 1e8;
bgdPhotonRate = 0;
dynamicRange = molPhotonRate - bgdPhotonRate;
pa = 0.166^2 ;
padsize=15;
md=3;
height=50;                  %length side of original image
block_width=2;                    %downsizing factor              

 
rho=pa*md/(block_width^2);  

%Generate Image
x  = molPhotonRate*generate_image( rho , height^2 );
gt_density=(block_width^2)*nnz(x)/(height^2*pa);

x=reshape(padarray(reshape(x,height,height),[padsize,padsize]),[],1);

padheight=height+2*padsize;
down_height=padheight/block_width;
n=padheight^2;


%Generate Super-resolution Matrix F
F = proj_sr_p( n, block_width,pa ,fwhm);



% Transformed image y
y=random('poiss',F*x)/dynamicRange;

x=x/dynamicRange;

delta = diag(max((F*x)/molPhotonRate,0.001/molPhotonRate));





%Options

%Options
opts.x0=x;
opts.t_max=30;
opts.eps_conv=1e-13;
opts.damp_meas = 0.5;
opts.height=padheight;
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


if mean((x_hat1 - x).^2)==0
    
    mse=0;
    recdensity=(block_width^2)*nnz(x_hat1)/(height^2*pa);
    fprintf('o The final MSE is %.4g\n', mse)   
    fprintf('o The recovered density is %.4g\n', recdensity)    
    fprintf('o The Ground truth density is %.4g\n', gt_density)    
    
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
    fprintf('o The final MSE is %.4g\n', mse)   
    fprintf('o The recovered density is %.4g\n', recdensity)    
    fprintf('o The Ground truth density is %.4g\n', gt_density)       
end