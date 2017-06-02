clear
seed=int64(14753);
addpath(genpath('..'));

molphotonrange=3:1:6;

rho_min=0.1; 
rho_step=0.3;
rho_max=11.8;  
md_vect=rho_min:rho_step:rho_max;  
mse_storm=zeros(length(md_vect),length(molphotonrange));
recdensity_storm=zeros(length(md_vect),length(molphotonrange));
mse_vamp=zeros(length(md_vect),length(molphotonrange));
recdensity_vamp=zeros(length(md_vect),length(molphotonrange));
gt_density=zeros(length(md_vect),length(molphotonrange));

t=1;

for j=molphotonrange
    
    
    
    
    
        molPhotonRate = 10^j;
        bgdPhotonRate = 0;
        dynamicRange = molPhotonRate - bgdPhotonRate;
        pa = 0.166^2 ;
        padsize=8;
        fwhm = 0.10;
   
        
        s=1; 
        height=64;                  %length side of original image
        block_width=2; 
        padheight=height+2*padsize;
        down_height=padheight/block_width;
        n=padheight^2;


        %Generate Super-resolution Matrix F
        F = proj_sr_p( n, block_width,pa ,fwhm);

        for md=md_vect

                                         %downsizing factor              


                        rho=pa*md/(block_width^2);  

                        %Generate Image
                        x  = molPhotonRate*generate_image( rho , height^2 );
                        gt_density(s,t)=(block_width^2)*nnz(x)/(height^2*pa);

                        x=reshape(padarray(reshape(x,height,height),[padsize,padsize]),[],1);


                        % Transformed image y
                        y=random('poiss',F*x);
                        figure(42);subplot(221);ima1 = imagesc(reshape(x, [padheight, padheight]));
                        colormap gray;axis image;pause(0.1);colorbar();title('original');
                        figure(42);subplot(222);ima2 = imagesc(reshape(y, [down_height, down_height]));
                        colormap gray;axis image;pause(0.1);colorbar();title('transform');
                        [ mse_storm(s,t), recdensity_storm(s,t) ] = fasterstorm( y,F,block_width,padheight,x,pa,height);
                        [ mse_vamp(s,t), recdensity_vamp(s,t) ] = vampoisson(y,F,block_width,padheight,x,pa,height,dynamicRange,n,down_height,rho);
               s=s+1;            
        end   
        
  

            figure(j); 
            plot(gt_density(:,t),mse_storm(:,t),'bo','DisplayName','Faster Storm')
            grid on;
            hold on;
            plot(gt_density(:,t),mse_vamp(:,t),'g*','DisplayName','Vamp')
            grid on;
            hold on;

        figure(j); 
        legend('Show','Location','SouthEast')
        xlabel('density')
        ylabel('Localization error (nm)')
         title(['Molecule Photon Rate = 1e',num2str(j)])
        hold off


         figure(j+length(molphotonrange)); 
            plot(gt_density(:,t),recdensity_storm(:,t),'bo','DisplayName','Faster Storm')
            grid on;
            hold on;
            plot(gt_density(:,t),recdensity_vamp(:,t),'g*','DisplayName','Vamp')
            grid on;
            hold on;
            plot(gt_density(:,t),gt_density(:,t),'--','DisplayName','Optimal')
        
        figure(j+length(molphotonrange)); 
        legend('Show','Location','SouthEast')
        xlabel('True Density')
        ylabel('Identified Density')
        title(['Molecule Photon Rate = 1e',num2str(j)])
        hold off
        
       t=t+1;
end


save('final_comparison_vamp_Faster_storm')