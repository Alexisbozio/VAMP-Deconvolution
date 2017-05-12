 function[Qx2,x2]=channel_multi(Ax2, Bx2, C, Az2, Bz2)
  
    A0=sparse(C')*sparse(diag(Az2))*sparse(C);
    B0=C'*Bz2;
    mat=A0+sparse(diag(Ax2));
    rhs = sparse(Bx2 +B0);
    
    tic
    Qx2=inv(mat);
    x2=(mat)\rhs;
    toc;
    
    