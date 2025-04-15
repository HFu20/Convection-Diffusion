function [ u_nPlusHf] = TVB_fun_Mstep_convec( u_nPlusHf, dt, matrix_A6, M, df, dx,nx,epsilon  )

mt = dt/M ;
lam = mt/(dx) ;
ubar = matrix_A6*u_nPlusHf ;
for i =1 : M

    %     t_k = T0+(k-1)*dt + (i-1)*mt ;
    %     t_kplus = T0+(k-1)*dt + (i)*mt ;

    
      [  fhat_m_plus_1,    fhat_m_minus_1 ] =  for_TVB_fhat_m( df, u_nPlusHf,dx, ubar,nx,epsilon ) ;
% u0 = 0 ; unx =0;
% [  fhat_m_plus_1,    fhat_m_minus_1 ] =  for_TVB_fhat_m_Diri( df, u_nPlusHf,dx, ubar,nx,epsilon,u0,unx ) ;
    u1 =  ubar - lam*( fhat_m_plus_1 - fhat_m_minus_1 )    ; %+mt*kxi ; %对流
    u_nPlusHfStar = Thomas(matrix_A6,u1 ) ;

    
     [  fhat_m_plus_2,    fhat_m_minus_2 ] =  for_TVB_fhat_m( df, u_nPlusHfStar,dx, ubar,nx,epsilon ) ;
%[  fhat_m_plus_2,    fhat_m_minus_2 ] =for_TVB_fhat_m_Diri( df, u_nPlusHfStar,dx, ubar,nx,epsilon,u0,unx ) ;
    ubar  =   0.5* ubar  + 0.5*u1 -   0.5*lam* (  fhat_m_plus_2 - fhat_m_minus_2)  ; % +mt*kxi ;   %对流

    u_nPlusHf =  Thomas(matrix_A6,ubar ) ;

end


