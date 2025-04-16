function [  u_nPlusHf4 ] = fun_Mstep_convec( u_nPlusHf2, matrix_A2x, matrix_Delta_x, lam_x, lam_y, matrix_Delta_y, matrix_A2y, dt,X, Y, gamma, T0,  k, kxi ,M )
 
for i =1 : M
    mt = dt/M ;
    t_k = T0+(k-1)*dt + (i-1)*mt ;
    t_kplus = T0+(k-1)*dt + (i)*mt ;
u_nPlusHf3    = u_nPlusHf2 - matrix_A2x\matrix_Delta_x*lam_x* fun_f(u_nPlusHf2)...
        - lam_y*(fun_g(u_nPlusHf2)) *matrix_Delta_y'/matrix_A2y' ...
        + mt*fun_Source (X,Y, gamma, t_k   )+ mt*kxi  ;  %g

    %-----------limiter------------------------------------------
    %    for i=1:nx
    %       u_nPlusHf3(i,:) = limiter_2sides_final(ny,u_nPlusHf3(i,:)',min(min(U0)),max(max(U0)))';
    %     end
    %---------------------------------------------------------------

    u_nPlusHf4    = u_nPlusHf2 -  matrix_A2x\matrix_Delta_x*lam_x* 0.5*(fun_f(u_nPlusHf2)+fun_f(u_nPlusHf3) ) ...%f
        -   lam_y/2* (fun_g(u_nPlusHf3)+fun_g(u_nPlusHf2) ) * matrix_Delta_y' /matrix_A2y' ...
        + mt/2 * ( fun_Source ( X,Y, gamma, t_k   ) + fun_Source ( X,Y, gamma,  t_kplus ) ) + mt*kxi ;%g;

end
end

