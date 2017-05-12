 function[a,c]=channel_sel(A, B, A0, B0)
       
      ll=1:length(A);
      %      tic
      mat=sparse(A0) + sparse(ll,ll,A);  
      rhs=sparse(B+B0);
       
      [c a]=sel2( mat, rhs ); 
     
    