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
height=60;                  %length side of original image
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
y=random('poiss',F*x);

[ mse, recdensity ] = fasterstorm( y,F,block_width,padheight,x,pa,height);

