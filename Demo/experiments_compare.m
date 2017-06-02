seed=int64(14753);
addpath(genpath('..'));
fwhm = 0.10;
molPhotonRate = 1e6;
bgdPhotonRate = 0;
dynamicRange = molPhotonRate - bgdPhotonRate;
pa = 0.166^2 ;
padsize=8;
height=64;                  %length side of original image
block_width=2;                    %downsizing factor  
s=20;    
rho_min=0.1; 
rho_step=0.3;
rho_max=12.7;  



md_vect=rho_min:rho_step:rho_max;    
mse_storm=zeros(size(md_vect));
recdensity_storm=zeros(size(md_vect));
mse_vamp=zeros(size(md_vect));
recdensity_vamp=zeros(size(md_vect));
gt_density=zeros(size(md_vect));
    


for md=5.8:rho_step:rho_max 
                
                           


                rho=pa*md/(block_width^2);  

                %Generate Image
                x  = molPhotonRate*generate_image( rho , height^2 );
                gt_density(s)=(block_width^2)*nnz(x)/(height^2*pa);

                x=reshape(padarray(reshape(x,height,height),[padsize,padsize]),[],1);

                padheight=height+2*padsize;
                down_height=padheight/block_width;
                n=padheight^2;


                %Generate Super-resolution Matrix F
                F = proj_sr_p( n, block_width,pa ,fwhm);

                % Transformed image y
                y=random('poiss',F*x);
                figure(42);subplot(221);ima1 = imagesc(reshape(x, [padheight, padheight]));
                colormap gray;axis image;pause(0.1);colorbar();title('original');
                figure(42);subplot(222);ima2 = imagesc(reshape(y, [down_height, down_height]));
                colormap gray;axis image;pause(0.1);colorbar();title('transform');
                [ mse_storm(s), recdensity_storm(s) ] = fasterstorm( y,F,block_width,padheight,x,pa,height);
                [ mse_vamp(s), recdensity_vamp(s) ] = vampoisson(y,F,block_width,padheight,x,pa,height,dynamicRange,n,down_height,rho);
       s=s+1;            
end   

    figure(1); 
    plot(md_vect,mse_storm,'DisplayName','Faster Storm')
    grid on;
    hold on;
    plot(md_vect,mse_vamp,'DisplayName','Vamp')
    grid on;
    hold on;
    
figure(1); 
legend('Show','Location','SouthEast')
xlabel('density')
ylabel('Localization error (nm)')
hold off


 figure(2); 
    plot(gt_density,recdensity_storm,'bo','DisplayName','Faster Storm')
    grid on;
    hold on;
    plot(gt_density,recdensity_vamp,'g*','DisplayName','Vamp')
    grid on;
    hold on;
    plot(gt_density,gt_density,'--','DisplayName','Optimal')
    figure(2); 
legend('Show','Location','SouthEast')
xlabel('True Density')
ylabel('Identified Density')
hold off
