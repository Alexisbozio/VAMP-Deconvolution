
height=60;


rho_min=0.02;
rho_max=0.2;
rho_step=0.02;
rho_size=ceil(1+(rho_max-rho_min)/rho_step);


airy_max=2;
airy_min=0.6;
step=0.2;
airy_size=ceil(1+(airy_max-airy_min)/step);


mse=zeros(airy_size,rho_size);
eltime=zeros(airy_size,rho_size);
t_final=zeros(airy_size,rho_size);

n=1;
for rho = rho_min:rho_step:rho_max
    
    g=num2str(rho);

   
 
    m=1;
    sigma_size=airy_min:step:airy_max;
    for airy = airy_min:step:airy_max
        
        fprintf('PSF=%.4g\n',airy)
        [mse(m,n) eltime(m,n) t_final(m,n)]=estimation_psf(rho,airy,height);
        
        m=m+1;
    end


    figure(1); 
    plot(sigma_size,eltime(:,n),'DisplayName',['rho=' num2str(rho)])
    grid on;
    hold on;
    
    figure(2); 
    plot(sigma_size,t_final(:,n),'DisplayName',['rho=' num2str(rho)])
    grid on;
    hold on;
    
    n=n+1;
    
end

figure(1); 
legend('Show','Location','SouthEast')
xlabel('psf size')
ylabel('elapse time')
hold off

figure(2); 
legend('Show','Location','SouthEast')
xlabel('psf size')
ylabel('number of iterations')
hold off


average_time_per_iteration=mean((eltime./t_final)');