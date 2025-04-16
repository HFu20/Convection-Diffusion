
%求空间阶数
clear
%------------------
N = [ 512 512 512 512];
N_t = [ 20 40 80 160 ] ;
% N = [   200 250 300 350 400];
gamma = 1/100;
df = 1 ;
cx = 1 ;  cy = 1 ;
  niu = 5/( 6*gamma ) ;

xmin = -1 ; xmax = 2 ;
ymin = -1 ; ymax = 2;

T0 = 0 ; T   = 0.5 ;
%%
nstep_x = length( N ) ;
Error_Table = zeros( nstep_x ,7 );
Bound_Table =  zeros( nstep_x ,3 );

for i = 1: nstep_x
    
    nx = N(i);
    ny = N(i);
  
dx = ( xmax - xmin )/ nx ;
dy = ( ymax - ymin ) / ny ;

X = xmin : dx : xmin + nx*dx ;
Y = ymin : dy : ymin + ny*dy ;
x = X( 1: end-1 ) ; % periodic
y = Y( 1: end-1 ) ;
% [ X1, Y1 ] = meshgrid( Y , X ) ; % (nx+1) * (nx+1)
[ X2, Y2 ] = meshgrid( y, x ) ; %  nx *  nx

% dt = dx/( 6*df ) ;
 nt = N_t( i );
% dt = dx*dx;
dt = ( ( T - T0 )/ nt ) ;
 T = T0 + nt*dt ;


%-- Initial ------------
U0 = Exact_sol( X2, Y2, cx, cy, T0, gamma ) ;
U  =  Exact_sol( X2, Y2, cx, cy, T, gamma ) ;
% u_n = U0  ;

%% common part
miu_x=dt/(dx*dx);
miu_y=dt/(dy*dy);
lam_x=dt/dx;
lam_y=dt/dy;

%  [E2,Emax,m_err,M_err] = fun_non2(dx,dt,xmin,xmax,gama,cx,T);
u = Solve_u_periodic( miu_x, miu_y,  lam_x, lam_y, nx , ny, gamma, nt, dt, X2, Y2, U0, T0 ) ; % Note X Y  
%-----------------------------------
% figure   
% mesh( X2, Y2, U ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("Exact solution")
% 
%  figure 
% mesh( X2, Y2, u ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("Numerical solution")
% 
%  figure 
% mesh( X2, Y2, U-u ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("Error")
%-----------------------
Emax = max(max(abs(U-u))) ;
E2 = sqrt(sum(dy*(sum(dx*((U-u).^2))))) ;
 E1 =  (sum(dy*(sum(dx*((U-u)))))) ;
M_err = max(max(U0))-max(max(u)) ;  %两个结果均为正，则说明 用了限制器 结果会保界
m_err = min(min(u))-min(min(U0)) ;


    % put n in the 1st column of table
    Error_Table(i,  1) = dt;
    % compute the L1 norm of error  and put it in the 4th column
    Error_Table(i, 2)  =  E2;
    % compute the L2 norm of error and put it in the 2nd column
    Error_Table(i, 4)  =  E2 ;
    % compute the Lmax norm of error and put it in the 2nd column
    Error_Table(i, 6)  =  Emax ;

    Bound_Table(i,  1) = dt;
    Bound_Table(i, 2) =  M_err;
    Bound_Table( i, 3 ) =  m_err;

    if i ~=1
        Error_Table(i, 3) = log(Error_Table(i, 2)/Error_Table(i-1, 2))...
            /log(Error_Table(i, 1)/Error_Table(i-1, 1));

        % compute the ratio for E2, and put it in the 5th column of table
        Error_Table(i, 5) = log(Error_Table(i, 4)/Error_Table(i-1, 4))...
            /log(Error_Table(i,1)/Error_Table(i-1, 1));

        % compute the ratio for Emax, and put it in the 5th column of table
        Error_Table(i, 7) = log(Error_Table(i, 6)/Error_Table(i-1, 6))...
            /log(Error_Table(i,1)/Error_Table(i-1, 1));
    end


end

format short e
% convert array to table for visualization
TableResults_Err = array2table(Error_Table, 'VariableNames', ...
    {'dt', 'E_1', 'ratio E_1', 'E_2', 'ratio E_2','E_max','ratio E_max'})
%--------------------------------------------------
TableResults_Bound = array2table(Bound_Table, 'VariableNames', {'dt', 'M_err', 'm_err'})





