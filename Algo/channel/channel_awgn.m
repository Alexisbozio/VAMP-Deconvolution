%  channel_awgn(A, B, channel_prmts; n_avgs = 1, eps_diff = 1e-3)
 
 %  Channel function...complicated because VAMP. 
 
  
  function[a,c]=channel_awgn(A, B, A0, B0,F,channel_prmts)
       
      ll=1:length(A);
      %      tic
      mat=sparse(A0) + sparse(ll,ll,A);  
      rhs=B+B0;
       
      %      con_num= cond( full(mat) )
      if(0)
      s = svd(full(mat));
      s=sort(s);
      s1=s(1:5)
      s2=s(end-5: end)
      end
   

       
      [c a]=sel2( mat, rhs ); 
      
        