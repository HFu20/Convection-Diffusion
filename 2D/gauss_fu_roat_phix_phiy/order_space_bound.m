%求空间阶数
clear
%------------------
%N_x = [ 50 100 150 200 ];
  N_x = [ 200 250 300 350];
gamma = 1e-3 ;
 %%
nstep_x = length( N_x ) ;
Error_Table = zeros( nstep_x ,5 );
Bound_Table =  zeros( nstep_x ,3 );

for i = 1: nstep_x
    
    nx = N_x(i);
 
%    [E2, Emax, wtime(i), m_err ,M_err ] = fun_order_Solve_u_inter_nonLag_split(nx, gamma ) ;
  [E2, Emax, wtime(i), m_err ,M_err ] = fun_order_Solve_u_inter_Lag_split(nx, gamma ) ;

    % put n in the 1st column of table
    Error_Table(i,  1) = nx;
    % compute the L1 norm of error  and put it in the 4th column
    Error_Table(i, 2)  =  E2;
    % compute the L2 norm of error and put it in the 2nd column
    Error_Table(i, 4)  =  Emax ;


    Bound_Table(i,  1) =  nx;
    Bound_Table(i, 2) =  M_err;
    Bound_Table( i, 3 ) =  m_err;

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
    {'Nx', 'E_2', 'ratio E_2','E_max','ratio E_max'})
%--------------------------------------------------
  TableResults_Bound = array2table(Bound_Table, 'VariableNames', {'Nx', 'm_err', 'M_err'})

writetable( TableResults_Err,'order_dtdx2_Lag.xls' ) ;
  writetable( TableResults_Bound,'BP_m_M_Lag.xls' ) ;





