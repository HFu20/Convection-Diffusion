clear

%% 23.04.21 画不同时刻的max_u min_u
% ================================================
% Eq: u_t+au_x=bu_{xx}
%-----------------------------------------------------------------
%exact solution: 
% u(x,y,t)= exp(-pi*pi*(ax+ay)*t)*sin(pi*(x-cx*t))*sin(pi*(y-cy*t))
% %------------参数----------------
%%%%  加和不加限制器的上下界
gamma=0.001; miu = 5/(3*gamma);df=1;
xmin =0;xmax = 3;
T = [   1 2 3 4 ] ;
T0=1;
nx = 200;
nt  = 200 ;
dx= (xmax  - xmin)/nx;
X=( xmin:dx:xmin+dx*nx)';
x=X(2:end-1); % 由于是周期 所以只需要一个边界点即可
cx = 1 ;

%===================================

Nt = length( T ) ;
min_u = zeros( Nt,1 ) ;
max_u = zeros( Nt,1 ) ;

[U0] = Exact_exp(T0,  X,  gamma) ;
min_u( 1 ) = min( U0 ) ;
max_u( 1 ) = max( U0 ) ;

for k = 2:  Nt

   T_final = T( k );
    dt = (T_final-T0)/nt;  %dt = dx/(3*df);
 
[ Iter,E1,E2,Emax,m_err,M_err,u] = fun_eg3_exp_diff_homo( dx,dt, nx, nt, T0, T_final, X,gamma ) ;
max_u( k ) = max(u ) ;
min_u( k ) = min( u ) ;
  

end

plot( T, max_u,'r')
hold on 
plot( T, min_u,'g' )
legend( 'max_u',' min_u ' )
xlabel( 'Time' )


