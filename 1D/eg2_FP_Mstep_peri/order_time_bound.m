 %测时间阶
clear
 
gamma = 1; miu = 5/(3*gamma); df=1;
xmin =-2*pi; xmax = 2*pi ;
T0=0;T =0.5;
N_t =  [ 100 200 400 800 1600  ];
D_t =  [ (T-T0)/N_t(1), (T-T0)/N_t(2), (T-T0)/N_t(3), (T-T0)/N_t(4), (T-T0)/N_t(5)  ] ;
%===================================
nstep_t = length( N_t )-1 ;
Error_Table = zeros( nstep_t ,5);
Bound_Table =  zeros( nstep_t ,3 ); 
K =1  ;
[ m_err_1,M_err_1,u_1] = fun_FP_t( N_t(1), D_t(1), T0, T,gamma,K,xmin,xmax   ) ;
[ m_err_2,M_err_2,u_2] = fun_FP_t( N_t(2),  D_t(2), T0, T,gamma,K, xmin,xmax ) ;
[ m_err_3,M_err_3,u_3] = fun_FP_t( N_t(3),  D_t(3), T0, T, gamma,K, xmin,xmax ) ;
[ m_err_4,M_err_4,u_4] = fun_FP_t( N_t(4),  D_t(4), T0, T, gamma,K, xmin,xmax ) ;
[ m_err_5,M_err_5,u_5] = fun_FP_t( N_t(5),  D_t(5), T0, T, gamma,K, xmin,xmax ) ;
 

 
E2(1) = sqrt(sum( D_t(1) *( ( u_1( 1:1:end  ) -  u_2( 1: 1: end ) ).^2))); 
E2(2) = sqrt(sum( D_t(2) *( ( u_2( 1: 1: end ) - u_3( 1: 1: end ) ).^2))); 
E2(3) = sqrt(sum( D_t(3) *( (  u_3( 1: 1: end ) -  u_4( 1: 1: end ) ).^2))); 
 E2(4) = sqrt(sum( D_t(4) *( (  u_4(  1: 1: end ) - u_5( 1: 1: end ) ).^2))); 

Emax(1) = max(max(abs( u_1( 1:1:end) -  u_2( 1: 1: end ) )  ));
Emax(2) = max(max(abs( u_2( 1: 1: end ) - u_3( 1: 1: end ) ) ));
Emax(3) = max(max(abs ( u_3( 1: 1: end ) -  u_4( 1: 1: end )  )));
 Emax(4) = max(max(abs(  u_4(  1: 1: end ) - u_5( 1: 1: end )  )));

    % 改用一阶Lagrang乘子法后 保正的同时，时间二阶，且保正
for i = 1:  nstep_t
    % put n in the 1st column of table
    Error_Table(i,  1) = N_t(i);
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
TableResults_Err_t = array2table(Error_Table, 'VariableNames', ...
    {'Nt','E_max', 'ratio E_max','E_2','ratio E_2'})
%--------------------------------------------------
% TableResults_Bound = array2table(Bound_Table, 'VariableNames', {'dx', 'M_err', 'm_err'})





