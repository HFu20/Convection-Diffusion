function [ u_nPlusHf] = fun_Mstep_convec( u_nPlusHf, dt, matrix_Delta_x,matrix_A2, kxi,T0, cx, x,gamma ,M,k) 

for i =1 : M
    mt = dt/M ;
    t_k = T0+(k-1)*dt + (i-1)*mt ;
    t_kplus = T0+(k-1)*dt + (i)*mt ;
%     I_x = eye( length(x),1 ) ; matrix_A2*  I_x*

    u_nPlusHfStar = u_nPlusHf + matrix_A2\( mt * matrix_Delta_x  *  fun_f(u_nPlusHf,x)  );%+ mt*kxi ;  %对流

%    u_nPlusHf  = u_nPlusHf + matrix_A2\(  0.5 *mt * matrix_Delta_x  ...
%                                 * (  fun_f(u_nPlusHf,x)  +  fun_f(u_nPlusHfStar,x) )  )   +  ...
%                                  mt*fun_Source(x, t_kplus, gamma,cx );%+mt*kxi ;   %对流
 u_nPlusHf  = 0.5*(u_nPlusHf + u_nPlusHfStar) + matrix_A2\(  0.5 *mt * matrix_Delta_x  ...
                                *    fun_f(u_nPlusHfStar,x)  ) ;   %对流
 
end

end

