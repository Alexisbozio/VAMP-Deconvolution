%  channel_awgn(A, B, channel_prmts; n_avgs = 1, eps_diff = 1e-3)

 %  Channel function...complicated because VAMP. 

 
 function[a,c]=channel_wood(A,B,A0, B0, F,delta)
  
   
 [m n]=size (F);
 I = sparse(eye(m)*delta);
 k=inv(sparse(diag(A)));
 C=k-k*F'*inv(I+F*k*F')*F*k;
 a=C*(B0+B);
 c=diag(C);