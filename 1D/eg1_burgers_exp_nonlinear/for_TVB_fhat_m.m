function  [ fhat_m_plus,  fhat_m_minus ] =  for_TVB_fhat_m( df, u,dx, ubar,nx,epsilon )

alpha = df ;

fhat_p = zeros( nx-1, 1 ) ;
fhat_n = zeros( nx-1, 1 ) ;
fbar_p = zeros( nx-1, 1 ) ;
fbar_n = zeros( nx-1, 1 ) ;
dfhat_p = zeros( nx-1, 1 ) ;
 dfhat_n = zeros( nx-1, 1 ) ;
delta_fbar_p = zeros( nx-1, 1 ) ;
delta_fbar_n = zeros( nx-1, 1 ) ;
Z_p = zeros( nx-1, 3 ) ;
Z_n = zeros( nx-1, 3 ) ;
dfhat_m_p = zeros( nx-1, 1 ) ;
dfhat_m_n = zeros( nx-1, 1 ) ;
fhat_p_m = zeros( nx-1, 1 ) ;
 fhat_n_m= zeros( nx-1, 1 ) ;
fhat_pm_minus = zeros( nx-1, 1 ) ;
fhat_nm_minus = zeros( nx-1, 1 ) ;

for i = 2 : nx-3
    fhat_p(i ) = 0.25*( fun_f( u(i+1) ) + fun_f( u(i) ) ) + 0.25*alpha*( u(i+1) + u(i) ) ;
    fbar_p(i) = 0.5*( fun_f( ubar(i) ) + alpha*ubar(i) ) ;
    dfhat_p(i) = fhat_p(i) - fbar_p(i) ;

    fhat_n(i ) = 0.25*( fun_f( u(i+1) ) + fun_f( u(i) ) ) - 0.25*alpha*( u(i+1) + u(i) ) ;
    fbar_n(i+1) = 0.5*( fun_f( ubar(i+1) ) - alpha*ubar(i+1) ) ;
    dfhat_n(i) = fbar_n(i+1) - fhat_n(i)   ;

    fbar_p(i+1) =0.5*( fun_f( ubar(i+1) ) + alpha*ubar(i+1) ) ;
    delta_fbar_p(i) = fbar_p(i+1) - fbar_p(i) ;
    fbar_p(i-1) = 0.5*( fun_f( ubar(i-1) ) + alpha*ubar(i-1) ) ;
    delta_fbar_p(i-1) =  fbar_p(i) - fbar_p(i-1) ;

    fbar_n(i+1) = 0.5*( fun_f(ubar(i+1)) - alpha*ubar(i+1)) ;
    fbar_n(i) = 0.5*( fun_f(ubar(i)) - alpha*ubar(i)) ;
    fbar_n(i+2) = 0.5*( fun_f(ubar(i+2)) - alpha*ubar(i+2)) ;
    delta_fbar_n(i) = fbar_n(i+1) - fbar_n(i) ;
    delta_fbar_n(i+1) = fbar_n(i+2) - fbar_n(i+1) ;

    Z_p(i,:) =  [ dfhat_p(i), delta_fbar_p(i), delta_fbar_p(i-1) ] ;
    Z_n(i,:) =  [ dfhat_n(i), delta_fbar_n(i), delta_fbar_n(i+1) ] ;
    dfhat_m_p(i) = minmod2( Z_p(i,:), epsilon, dx ) ;
    dfhat_m_n(i) = minmod2( Z_n(i,:), epsilon, dx ) ;
    fhat_p_m(i) = fbar_p(i) +dfhat_m_p(i) ;
    fhat_n_m(i) = fbar_n(i+1) - dfhat_m_n(i) ;

end

i = 1 ;
fhat_p(i) = 0.25*( fun_f(u(i+1)) + fun_f(u(i)) + alpha*( u(i+1)+ u(i) )) ;
%fbar_p(i) = 0.5*( fun_f( ubar(i) ) + alpha*ubar(i) ) ;
fbar_p_0 = 0.5*( fun_f( 1/6*u(i) ) + alpha*1/6*u(1) ) ;
dfhat_p(i) = fhat_p(i) - fbar_p(i) ;
delta_fbar_p(i) = fbar_p(i+1) - fbar_p(i) ;
delta_fbar_p_0 = fbar_p(i)  -  fbar_p_0;
Z_p(i,:) =  [ dfhat_p(i), delta_fbar_p(i), delta_fbar_p_0 ] ;
dfhat_m_p(i) =  minmod2( Z_p(i,:), epsilon, dx ) ;
fhat_p_m(i) = fbar_p(i) + dfhat_m_p(i);

fhat_n(i)  = 0.25*( fun_f(u(i+1)) + fun_f(u(i)) - alpha*( u(i+1)+ u(i) ) ) ;
fbar_n(i) = 0.5*( fun_f(ubar(i)) - alpha*ubar(i) ) ;
 dfhat_n(i) = fbar_n(i+1) - fhat_n(i)   ; 
 delta_fbar_n(i) = fbar_n(i+1) - fbar_n(i) ;
 delta_fbar_n(i+1) = fbar_n(i+2) - fbar_n(i+1) ;
 Z_n(i,:) =  [ dfhat_n(i), delta_fbar_n(i), delta_fbar_n(i+1) ] ;
 dfhat_m_n(i) = minmod2( Z_n(i,:), epsilon, dx ) ;
 fhat_n_m(i) = fbar_n(i+1) - dfhat_m_n(i) ;


i = nx-2;

fbar_p(i+1) = 0.5*( fun_f(ubar(i+1)) + alpha*ubar(i+1) ) ;
fbar_n(i+1) = 0.5*( fun_f( ubar(i+1) ) - alpha*ubar(i+1) ) ;
fbar_n(i+2) = 0.5*( fun_f( 1/6*ubar(i+1) ) -  alpha*1/6*ubar(i+1) ) ;

fhat_p(i) = 1/4*(fun_f(u(i)) +fun_f(u(i+1))) + alpha/4*( u(i) + u(i+1)) ;
fhat_n(i) =  1/4*(fun_f(u(i)) +fun_f(u(i+1))) - alpha/4*( u(i) + u(i+1)) ;

dfhat_p(i) = fhat_p(i) - fbar_p(i) ;
dfhat_n(i) = fbar_n(i+1) - fhat_n(i) ;
delta_fbar_p(i) = fbar_p(i+1) - fbar_p(i) ;
delta_fbar_n(i) = fbar_n(i+1) - fbar_n(i) ;
delta_fbar_n(i+1) = fbar_n( i+1 ) - fbar_n(i) ;

Z_p(i,:) = [ dfhat_p(i), delta_fbar_p(i), delta_fbar_p(i-1) ] ;
Z_n( i,: ) = [ dfhat_n(i), delta_fbar_n(i), delta_fbar_n(i+1) ] ;
dfhat_m_p(i) =  minmod2( Z_p(i,:), epsilon, dx ) ;
fhat_p_m(i) = fbar_p(i) + dfhat_m_p(i);
dfhat_m_n(i) =  minmod2( Z_n(i,:), epsilon, dx ) ;
fhat_n_m(i) = fbar_n(i+1) - dfhat_m_n(i);


i = nx-1 ;
fhat_p(i) = 0.25*( fun_f(u(i)) + alpha*u(i) ) ;
fbar_p(i+1) = 0.5*( fun_f(1/6*u(i)) + alpha*1/6*u(i) ) ;
dfhat_p(i) = fhat_p(i) - fbar_p(i) ;
fhat_n(i) = 1/4*( fun_f(u(end)) - alpha*u(end) ) ;
fbar_n(i+1) = 1/2*( fun_f(1/6*u(end)) -alpha/6*u(end) ) ;

delta_fbar_p(i) = fbar_p(i+1) - fbar_p(i) ;
delta_fbar_n(i) = fbar_n(i+1) - fbar_n(i) ;
delta_fbar_n(i+1) = -fbar_n( i+1 );
dfhat_p(i) = fhat_p(i) - fbar_p(i) ;
dfhat_n(i) = fbar_n(i+1) - fhat_n(i) ;

Z_p(i,:) = [ dfhat_p(i), delta_fbar_p(i), delta_fbar_p(i-1) ] ;
Z_n( i,: ) = [ dfhat_n(i), delta_fbar_n(i), delta_fbar_n(i+1) ] ;
dfhat_m_p(i) =  minmod2( Z_p(i,:), epsilon, dx ) ;
fhat_p_m(i) = fbar_p(i) + dfhat_m_p(i);
dfhat_m_n(i) =  minmod2( Z_n(i,:), epsilon, dx ) ;
fhat_n_m(i) = fbar_n(i+1) - dfhat_m_n(i); 

 

%--------------------------
 fhat_pm_plus =  fhat_p_m ;
 fhat_pm_minus( 2:end ) =  fhat_p_m( 1:end-1 ) ;
  fhat_nm_plus =  fhat_n_m ;
  fhat_nm_minus( 2:end ) =  fhat_n_m( 1:end-1 ) ;

  %=========================
  fhat_p_0 = 0.25*( fun_f(u(1)) + alpha*u(1) ) ;
  %fbar_p_0 = 0.5*( fun_f(1/6*u(1)) + alpha*1/6*u(1) ) ;
  dfhat_p_0 = fhat_p_0 - fbar_p_0 ;
  delta_fbar_p_0 = fbar_p(1) -  fbar_p_0 ;
   delta_fbar_p_minus =    fbar_p_0 ;
   Z_p_0 =  [ dfhat_p_0, delta_fbar_p_0, delta_fbar_p_minus ] ;
dfhat_m_p_0 =  minmod2( Z_p_0, epsilon, dx ) ;
fhat_p_m_0 = fbar_p_0 + dfhat_m_p_0 ; 

  fhat_n_0 = 0.25*( fun_f(u(1)) - alpha*u(1) ) ;
  fbar_n_0 = 0.5*( fun_f(1/6*u(1)) -alpha/6*u(1) ) ;
  dfhat_n_0 = fbar_n(1) - fhat_n_0 ;
  delta_fbar_n_0 = fbar_n(1) -  fbar_n_0 ;
  %delta_fbar_n
  Z_n_0 =  [ dfhat_n_0, delta_fbar_n_0, delta_fbar_n(1) ] ;
 dfhat_n_m_0  = minmod2( Z_n_0, epsilon, dx ) ;
fhat_n_m_0 = fbar_n(1) - dfhat_n_m_0 ; 
 
  fhat_pm_minus(1) = fhat_p_m_0  ;
fhat_nm_minus(1) = fhat_n_m_0;

   fhat_m_plus =  fhat_pm_plus + fhat_nm_plus ;
 fhat_m_minus =  fhat_pm_minus + fhat_nm_minus ;
 


end