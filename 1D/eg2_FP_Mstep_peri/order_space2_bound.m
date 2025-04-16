%求空间阶数
% 取 dt = miu* dx^2 同时验证了时间精度
%% 23.04.08
% 四阶空间精度
clear

%%================================================
% Eq: u_t+au_x=bu_{xx}
%-----------------------------------------------------------------
%exact solution: 
% u(x,y,t)= exp(-pi*pi*(ax+ay)*t)*sin(pi*(x-cx*t))*sin(pi*(y-cy*t))
% %------------参数----------------
%%%%  加和不加限制器的上下界
gamma = 1; miu = 5/(3*gamma); df=1;
xmin =-2*pi; xmax = 2*pi ;
T0=0;T =0.2;
N_x = [ 200 400 800 1600 3200 6400];
% D_x =  [ (xmax-xmin)/N_x(1), (xmax-xmin)/N_x(2), (xmax-xmin)/N_x(3), (xmax-xmin)/N_x(4), (xmax-xmin)/N_x(5) ] ;
% D_x =  [ (xmax-xmin)/N_x(1), (xmax-xmin)/N_x(2), (xmax-xmin)/N_x(3), (xmax-xmin)/N_x(4) ] ;
D_x =  [ (xmax-xmin)/N_x(1), (xmax-xmin)/N_x(2), (xmax-xmin)/N_x(3), (xmax-xmin)/N_x(4), (xmax-xmin)/N_x(5) ,(xmax-xmin)/N_x(6) ] ;
% N_t = [  750 850 950 1050 1150 ];
x1 = (xmin:D_x(1):xmax)' ; 
x2 = (xmin:D_x(2):xmax)' ; 
x3 = (xmin:D_x(3):xmax )';
x4 = (xmin:D_x(4):xmax)' ;
x5 = (xmin:D_x(5):xmax)' ;
x6 = (xmin:D_x(6):xmax)' ;
%===================================

nstep_x = length( N_x )-1 ;
Error_Table = zeros( nstep_x ,5);
Bound_Table =  zeros( nstep_x ,3 );
% K =10 ;
K =1 ;
[ m_err_1,M_err_1,u_1] = fun_FP( N_x(1), D_x(1), T0, T,gamma,K,xmin,xmax   ) ;
[ m_err_2,M_err_2,u_2] = fun_FP( N_x(2),  D_x(2), T0, T,gamma,K, xmin,xmax ) ;
[ m_err_3,M_err_3,u_3] = fun_FP( N_x(3),  D_x(3), T0, T, gamma,K, xmin,xmax ) ;
[ m_err_4,M_err_4,u_4] = fun_FP( N_x(4),  D_x(4), T0, T, gamma,K, xmin,xmax ) ;
[ m_err_5,M_err_5,u_5] = fun_FP( N_x(5),  D_x(5), T0, T, gamma,K, xmin,xmax ) ;
[ m_err_6,M_err_6,u_6] = fun_FP( N_x(6),  D_x(6), T0, T, gamma,K, xmin,xmax ) ;


E2(1) = sqrt(sum( D_x(1) *( (   u_1( 1:1:end ) -  u_6( 1: 32: end ) ).^2))); 
E2(2) = sqrt(sum( D_x(2) *( (  u_2( 1: 2: end ) - u_6( 1: 32: end ) ).^2))); 
E2(3) = sqrt(sum( D_x(3) *( (  u_3( 1: 4: end ) -  u_6( 1: 32: end ) ).^2))); 
 E2(4) = sqrt(sum( D_x(4) *( (  u_4(  1: 8: end ) - u_6( 1: 32: end ) ).^2))); 
 E2(5) = sqrt(sum( D_x(4) *( (  u_5(  1: 16: end ) - u_6( 1: 32: end ) ).^2))); 

Emax(1) = max(max(abs(  u_1( 1:1:end ) -  u_6( 1: 32: end ) )  ));
Emax(2) = max(max(abs(  u_2( 1: 2: end ) - u_6( 1: 32: end ) ) ));
Emax(3) = max(max(abs(  u_3( 1: 4: end ) -  u_6( 1: 32: end )  )));
Emax(4) = max(max(abs(  u_4(  1: 8: end ) - u_6( 1: 32: end )  )));
Emax(5) = max(max(abs(  u_5(  1: 16: end ) - u_6( 1: 32: end )  )));
% [E1, E2,Emax,m_err,M_err] = fun_t_non_Lag_sin(T0,T,dx,dt,nx,nt, x,gamma,cx);

    % 改用一阶Lagrang乘子法后 保正的同时，时间二阶，且保正
for i = 1:  nstep_x
    % put n in the 1st column of table
    Error_Table(i,  1) = N_x(i);
    % compute the L1 norm of error  and put it in the 4th column
    Error_Table(i, 2)  =  Emax(i) ;
    % compute the Lmax norm of error and put it in the 2nd column
    Error_Table(i, 4)  =  E2(i) ;

%     Bound_Table(i, 2) =  M_err_2;
%     Bound_Table( i, 3 ) =  m_err_2;

    if i ~=1
        Error_Table(i, 3) =- log(Error_Table(i, 2)/Error_Table(i-1, 2))...
            /log(Error_Table(i, 1)/Error_Table(i-1, 1));

        % compute the ratio for E2, and put it in the 5th column of table
        Error_Table(i, 5) =- log(Error_Table(i, 4)/Error_Table(i-1, 4))...
            /log(Error_Table(i,1)/Error_Table(i-1, 1));

    end


end

format short e
% convert array to table for visualization
TableResults_Err = array2table(Error_Table, 'VariableNames', ...
    {'Nx','E_max', 'ratio E_max','E_2','ratio E_2'})
%--------------------------------------------------
% TableResults_Bound = array2table(Bound_Table, 'VariableNames', {'dx', 'M_err', 'm_err'})





