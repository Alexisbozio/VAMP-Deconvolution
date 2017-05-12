function [z sol] = selinv(A, b) 
   A=tril(A);
   [z sol] = florent(A,b);
