function[a1, A1, B1]= vamp_stream_wood(y, F,Ap,Bp, opts)
    prior = opts.prior;
    prior_prmts = opts.prior_prmts;
    channel_prmts = opts.channel_prmts;
    channel=opts.channel;
    x0 = opts.x0;               %ground truth
    t_max = opts.t_max;
    eps_conv = opts.eps_conv;
    damp_meas = opts.damp_meas;
  
     
    [m n]=size(F);

    
  
    
    A1 = zeros(n,1);
    B1 = zeros(n,1);
    A2 = zeros(n,1);
    B2 = zeros(n,1);
    a1 = zeros(n,1);
    c1 = ones(n,1);
    
    F=sparse(F);
    y=sparse(y);
    channel_prmts=sparse(channel_prmts);
   
    A0 = F' * F ./ channel_prmts;
    B0 = F' * y ./ channel_prmts;


    for t = 1:t_max
        a1_old = a1;

        %A2_new = 1. ./ c1 - A1
        A2_new = max(1. ./ c1 - A1, 1e-11);
        B2_new = a1 ./ c1 - B1;
        A2 = damp(A2_new, A2,damp_meas); 
        B2 = damp(B2_new, B2,damp_meas);
        [a2, c2] = channel(A2,B2,A0, B0, F,channel_prmts);
        
               
        %A1_new = 1. ./ c2 - A2
        A1_new = max(1. ./ c2 - A2, 1e-11);
        B1_new = a2 ./ c2 - B2;
        A1 = damp(A1_new, A1,damp_meas); 
        B1 = damp(B1_new, B1,damp_meas);
        [a1, c1] = prior(A1+Ap, B1+Bp, prior_prmts);

           
       
      

        diff = mean((a1 - a1_old).^2);
        rss = mean((y - F * a1).^2);
        %if x0 ~= nothing 
            mse = mean((a1 - x0).^2); 
        %end
        fprintf( ' t=%.4g, diff = %.4g; rss = %.4g, mse = %.4g\n', t, diff, rss, mse)
        if diff < eps_conv 
            break 
        end
        
        
         
            
       
    end
    
    
    
