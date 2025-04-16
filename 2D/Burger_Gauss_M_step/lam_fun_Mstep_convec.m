function [  u_nPlusHf4 ] = lam_fun_Mstep_convec( u_nPlusHf2, matrix_A2x, matrix_Delta_x, lam_x, lam_y, matrix_Delta_y, matrix_A2y, dt,X, Y, gamma, T0,  k,  M)
 
% lam_fun_Mstep_convec( u_nPlusHf2, matrix_A2x, matrix_Delta_x, lam_x, lam_y, matrix_Delta_y, matrix_A2y, dt,X, Y, gamma, T0,  k,  M,kxi  ,lam_n, u_n ,gprime )
 
for i =1 : M
    mt = dt/M ;
    t_k = T0+(k-1)*dt + (i-1)*mt ;
    t_kplus = T0+(k-1)*dt + (i)*mt ;
     u_nPlusHf3    = u_nPlusHf2 - matrix_A2x\matrix_Delta_x*lam_x* fun_f(u_nPlusHf2)...
        - lam_y*(fun_g(u_nPlusHf2)) *matrix_Delta_y'/matrix_A2y' ...
        + mt*fun_Source (X,Y, gamma, t_k   )  ; %  +  mt*(kxi + lam_n.*gprime(u_n));   %g

    u_nPlusHf4    = u_nPlusHf2 -  matrix_A2x\matrix_Delta_x*lam_x* 0.5*(fun_f(u_nPlusHf2)+fun_f(u_nPlusHf3) ) ...%f
        -   lam_y/2* (fun_g(u_nPlusHf3)+fun_g(u_nPlusHf2) ) * matrix_Delta_y' /matrix_A2y' ...
        + mt/2 * ( fun_Source ( X,Y, gamma, t_k   ) + fun_Source ( X,Y, gamma,  t_kplus ) )  ; % + 0.5* mt*(kxi + lam_n.*gprime(u_n)); %g;
%% ADI
%   un3_hat    = u_nPlusHf2*matrix_A2y' - matrix_A2x\matrix_Delta_x*lam_x* fun_f(u_nPlusHf2)*matrix_A2y'...
%                         - lam_y*(fun_g(u_nPlusHf2)) *matrix_Delta_y' ...
%                         + mt*fun_Source (X,Y, gamma, t_k   ) *matrix_A2y' ;  %g
%   u_nPlusHf3 = un3_hat/matrix_A2y' ;
%   u_n4_hat    = u_nPlusHf2*matrix_A2y'  -  matrix_A2x\matrix_Delta_x*lam_x* 0.5*(fun_f(u_nPlusHf2)+fun_f(u_nPlusHf3) )*matrix_A2y'  ...%f
%         -   lam_y/2* (fun_g(u_nPlusHf3)+fun_g(u_nPlusHf2) ) * matrix_Delta_y'  ...
%         + mt/2 * ( fun_Source ( X,Y, gamma, t_k   ) + fun_Source ( X,Y, gamma,  t_kplus ) )  *matrix_A2y'; %+ mt*kxi ;%g;
%    u_nPlusHf4 = u_n4_hat/matrix_A2y' ;
end
end

