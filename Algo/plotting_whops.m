
for n=1:2
    rho=n/10;


    figure(1); 
    plot(sigma_size,eltime(:,n),'DisplayName',['rho=' num2str(rho)])
    grid on;
    hold on;
    
    figure(2); 
    plot(sigma_size,t_final(:,n),'DisplayName',['rho=' num2str(rho)])
    grid on;
    hold on;
end
%{
n=2;
rho=n/10;

figure(1); 
plot(sigma_size,eltime(:,n),'DisplayName',['rho=' num2str(rho)])
%}

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