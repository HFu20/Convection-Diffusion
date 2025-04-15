function [u]=Thomas(A,d)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  * @file   Thomas.m
%  * @author Fu hongfei
%  * @date    2016/02
%  * 
%  * @brief  Thomas algorithm for solving tridiagonal matrix system[a_i;b_i;c_i]
%  * input   a, lower diagonal; 
%  *         b, main diagonal; 
%  *         c, upper diagonal
%  *         d, right-hand side vector
%  *         N, size of vector
%  * output  u, the solution of Euler method at t+dt
%% Thomas alogrithm for TDMA system
 a = diag(A,-1) ;
 b = diag(A) ;
 c = diag(A,1) ;
    % Forward elimination phase
    N=length(d);
    u=zeros(N,1);
    for k=2:N
        m=a(k-1)/b(k-1);
        b(k)=b(k) - m*c(k-1);  
        d(k)=d(k) - m*d(k-1); 
    end
%% Backward substitution phase
    u(N)= d(N)/b(N);
    for k=N-1:-1:1 
        u(k)=(d(k)- c(k)* u(k+1))/b(k);  
    end
end