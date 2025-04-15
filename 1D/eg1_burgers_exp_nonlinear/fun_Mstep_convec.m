function [ u_nPlus1_tilde] = fun_Mstep_convec( u_nPlusHf, dt, matrix_Delta_x,matrix_A2, kxi,T0, cx,x,gamma ,M,k )

kxi = 0 ;
for i =1 : M
    mt = dt/M ;
    t_k = T0+(k-1)*dt + (i-1)*mt ;
    t_kplus = T0+(k-1)*dt + (i)*mt ;
    I_x = eye( length(x),1 ) ;

    u_nPlusHfStar = u_nPlusHf - matrix_A2\( mt * matrix_Delta_x  *  fun_f(u_nPlusHf)   ) + ...
                                   mt*fun_Source(x, t_k, gamma,cx ) ; %+mt*kxi ; %对流

    u_nPlus1_tilde  = u_nPlusHf - matrix_A2\(  0.5 *mt * matrix_Delta_x  ...
                                * (  fun_f(u_nPlusHf)  +  fun_f(u_nPlusHfStar) )  )   +  ...
                                 mt/2*( fun_Source(x, t_k, gamma,cx ) + fun_Source(x, t_kplus, gamma,cx ) ) ; % +mt*kxi ;   %对流
   u_nPlusHf =  u_nPlus1_tilde ;
end

end

 
