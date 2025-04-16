
%求空间阶数
clear
%------------------
% N_x = [ 50 100 150 200 ];
%  N_x = [   50 100 150 200];
% gamma = 1/64;
% N_x = [ 200, 250,300,350];
%  N_x = [ 150,200, 250,300];
N_x = [ 90,180,270,360];
gamma = 0.005;
df = 1 ;
cx = 1 ;  cy = 1 ;
niu = 5/( 6*gamma ) ;

xmin = -1 ; xmax = 2 ;
ymin = -1 ; ymax = 2;

T0 = 0 ; T_old   = 0.6 ;
%%
nstep_x = length( N_x ) ;
Error_Table_BDF = zeros( nstep_x ,5 );
Error_Table_CN = zeros( nstep_x ,5 );

for i = 1: nstep_x

    nx = N_x(i);
    ny = N_x(i);

    dx = ( xmax - xmin )/ nx ;
    dy = ( ymax - ymin ) / ny ;

    X = xmin : dx : xmin + nx*dx ;
    Y = ymin : dy : ymin + ny*dy ;
    x = X( 2: end-1 ) ;
    y = Y(2: end-1 ) ;
    [ X1, Y1 ] = meshgrid( Y , X ) ; % (nx+1) * (nx+1)
    [ X2, Y2 ] = meshgrid( y, x ) ; %

    % dt = dx/( 6*df ) ;
    dt =dx*dx;
    % dt = dx*dx;
    nt = floor( ( T_old - T0 )/ dt ) ;
    T = T0 + nt*dt ;
if T ~= T_old
    fprintf('T_new = %f',T)
end

    %-- Initial ------------
    U0 = Exact_sol( X1, Y1, cx, cy, T0, gamma ) ;
    U  =  Exact_sol( X1, Y1, cx, cy, T, gamma ) ;
    % u_n = U0  ;

    %% common part
    miu_x=dt/(dx*dx);
    miu_y=dt/(dy*dy);
    lam_x=dt/dx;
    lam_y=dt/dy;


    %-----------------------------------
    [u_BDF ] = fun_bdf2_2d( miu_x, miu_y,  lam_x, lam_y, nx , ny,  gamma, nt, dt, X2, Y2,U0, T0) ;
    %-----------------------
    Emax_BDF = max(max(abs(U-u_BDF))) ;
    E2_BDF = sqrt(sum(dy*(sum(dx*((U-u_BDF).^2))))) ;

    % put n in the 1st column of table
    Error_Table_BDF(i,  1) = dx;
    % compute the L2 norm of error and put it in the 2nd column
    Error_Table_BDF(i, 2)  =  Emax_BDF ;
    % compute the Lmax norm of error and put it in the 2nd column
    Error_Table_BDF(i, 4)  =  E2_BDF ;

   

    if i ~=1
        % compute the ratio for E2 and put it in the 5th column of table
        Error_Table_BDF(i, 3) = log(Error_Table_BDF(i, 2)/Error_Table_BDF(i-1, 2))...
            /log(Error_Table_BDF(i, 1)/Error_Table_BDF(i-1, 1));

        % compute the ratio for Emax and put it in the 5th column of table
        Error_Table_BDF(i, 5) = log(Error_Table_BDF(i, 4)/Error_Table_BDF(i-1, 4))...
            /log(Error_Table_BDF(i,1)/Error_Table_BDF(i-1, 1)); 
    end


end



format short e
% convert array to table for visualization
TableResults_Err_BDF = array2table(Error_Table_BDF, 'VariableNames', ...
    {'dx', 'E_max', 'ratio E_max','E_2','ratio E_2'})
 

