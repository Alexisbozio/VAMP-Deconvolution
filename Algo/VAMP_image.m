function [mse,a1] = VAMP_image(y,F,opts)

    channel = opts.channel;
    channel_prmts = opts.channel_prmts;
    x0 = opts.x0;              
    t_max = opts.t_max;
    eps_conv = opts.eps_conv;
    damp_meas = opts.damp_meas;
    height=opts.height;
    width=opts.width;
  
   [m n] = size(F);

 
    A1 = zeros(n,1);
    B1 = zeros(n,1);
    A2 = zeros(n,1);
    B2 = zeros(n,1);
    a1 = zeros(n,1);
    c1= ones(n,1);

    
    F=sparse(F);
    y=sparse(y);
    channel_prmts=sparse(channel_prmts);
   
    A0 = F' * F ./ channel_prmts;
    B0 = F' * y ./ channel_prmts;
    mse=zeros(1,t_max);

    for t = 1:t_max
        a1_old = a1;

        
        
        %A2_new = 1. ./ c1 - A1
        A2_new = max(1. ./ c1 - A1, 1e-11);
        B2_new = a1 ./ c1 - B1;
        A2 = damp(A2_new, A2,damp_meas); 
        B2 = damp(B2_new, B2,damp_meas);
        [a2, c2] = channel(A2, B2, A0,B0);
        c2=mean(c2);
        fprintf('    - mean(c2): %0.3f\n',c2);
        figure(42);
            subplot(221);
            imagesc(reshape(a2,height,height));
            colormap gray;
            axis image;
            title('a2');
            colorbar();
            pause(0.1);
        
        
        %A1_new = 1. ./ c2 - A2
        A1_new = max(1. ./ c2 - A2, 1e-11);
        B1_new = a2 ./ c2 - B2;
        A1 = damp(A1_new, A1,damp_meas); 
        B1 = damp(B1_new, B1,damp_meas);
        r1=255.*B1./A1;
        v1=mean(1./A1);
        
        R1=reshape(r1,[width,height]);
        [NA, a1]=BM3D(1,R1,v1,'np',0);
        a1=a1(:);
        fprintf('    - mean(v1): %0.3f\n',v1);
        
        figure(42);
            subplot(222);
            imagesc(reshape(a1,height,height));
            colormap gray;
            axis image;
            title('a1');
            colorbar();
            pause(0.1);
        
        
        eta=randn(n,1);
        epsilon=max(r1(:))/1000+eps;
        
        l1=r1+epsilon.*eta;
        l1=reshape(l1,[width,height]);
        [NA, g]=BM3D(1,l1,v1,'np',0);
        g=g(:);
        c1=eta'*(g-a1)/epsilon;
        fprintf('    - c1: %0.3f\n',c1);

        diff = mean((a1 - a1_old).^2);
        rss = mean((y - F * a1).^2);
        mse(t) = mean((a1 - x0).^2); 
     
        fprintf( ' t=%.4g, diff = %.4g; rss = %.4g, mse = %.4g\n', t, diff, rss, mse(t))
        if diff < eps_conv 
            break 
        end
     
       
    end
    

