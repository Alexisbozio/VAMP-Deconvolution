function[ax1]= vamp_corean(y, C,D, opts)
    prior = opts.prior;
    prior_prmts = opts.prior_prmts;
    channel_prmts = opts.channel_prmts;
    x0 = opts.x0;               %ground truth
    t_max = opts.t_max;
    eps_conv = opts.eps_conv;
    damp_meas = opts.damp_meas;
    height=opts.height;
    %plot = opts.plot;

    [m n]=size(D);

    %First Layer
    Ax1 = zeros(n,1);
    Bx1 = zeros(n,1);
    Ax2 = ones(n,1);
    Bx2 = zeros(n,1);
    
  
    %Second Layer
    Az1 = ones(n,1);
    Bz1 = zeros(n,1);
    Az2 = zeros(n,1);
    Bz2 = zeros(n,1);
    
    ax1= randn(n,1);
    A0 = sparse(D' * D ./ channel_prmts);
    B0 = sparse(D' * y ./ channel_prmts);


    for t = 1:t_max
        a1_old = ax1;
        
      
        %{
        %In z1 ; receive Az1 Bz1
        %Compute az1 cz1
        w = Bz1 ./ Az1;  
        v = 1 ./ Az1;
        [g, dg] = channel_probit(w, v, y,channel_prmts);
        
        cz1 = v + v.^2 .*dg;
        az1 = w + v.*g;
        %} 
        
        %A2_new = 1. ./ c1 - A1

       [az1, cz1] = channel_awgn(Az1, Bz1, A0,B0,C,channel_prmts);
        
               
        
        %Compute Az2 Bz2
        Az2_new = 1. ./ cz1 - Az1;
        Bz2_new = az1 ./ cz1 - Bz1;
        Az2 = damp(Az2_new, Az2,damp_meas); 
        Bz2 = damp(Bz2_new, Bz2,damp_meas);
        
        
        %In x2 ; receive Ax2 Bx2 & Az2 Bz2
        % Compute ax2 cx2
        [Qx2 x2]=channel_multi(Ax2,Bx2,C,Az2,Bz2);
        %Compute Ax1, Bx1
        cx2=diag(Qx2);
        Ax1_new = 1. ./ cx2 - Ax2;
        Bx1_new = x2 ./ cx2 - Bx2;
        Ax1 = damp(Ax1_new, Ax1,damp_meas); 
        Bx1 = damp(Bx1_new, Bx1,damp_meas);
        
        
        %In x1 ; receive Ax1 Bx1
        %Compute cx1, ax1
        [ax1, cx1] = prior(Ax1, Bx1,prior_prmts);
        %Compute Ax2, Bx2
        Ax2_new=1. ./ cx1 - Ax1;
        Bx2_new = ax1 ./ cx1 - Bx1;
        Ax2 = damp(Ax2_new, Ax2,damp_meas); 
        Bx2 = damp(Bx2_new, Bx2,damp_meas);
        
        
        
        figure(42);
            subplot(121);
            imagesc(reshape(ax1,height,height));
            colormap gray;
            axis image;
            title('ax1');
            colorbar();
            pause(0.1); 
        
            
        %In z2 ; receive Az2 Bz2 & ax2 cx2
        %Compute cz2 z2
        z2=C*x2;
        Qz2=sparse(C)*sparse(Qx2)*sparse(C)';
        
        cz2=diag(Qz2);
        %Compute Az1 Bz1
        Az1_new=1. ./ cz2 - Az2;
        Bz1_new = z2 ./ cz2 - Bz2;
        Az1 = damp(Az1_new, Az1,damp_meas); 
        Bz1 = damp(Bz1_new, Bz1,damp_meas);
        
        figure(42);
            subplot(122);
            imagesc(reshape(z2,height,height));
            colormap gray;
            axis image;
            title('z2');
            colorbar();
            pause(0.1);   
       
      

        diff = mean((ax1 - a1_old).^2);
        rss = mean((y - D*C * ax1).^2);
        %if x0 ~= nothing 
            mse = mean((ax1 - x0).^2); 
        %end
        fprintf( ' t=%.4g, diff = %.4g; rss = %.4g, mse = %.4g\n', t, diff, rss, mse)
        if diff < eps_conv 
            break 
        end

       
    end
    
    
    
    
