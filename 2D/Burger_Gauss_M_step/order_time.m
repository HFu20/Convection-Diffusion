
%求空间阶数
clear
%------------------
N = [ 512 512 512 512];
N_t = [ 20 40 80 160 ] ;
% N = [   200 250 300 350 400];
gamma = 1/64;
df = 1 ;
cx = 1 ;  cy = 1 ;
  niu = 5/( 6*gamma ) ;

xmin = -1 ; xmax = 2 ;
ymin = -1 ; ymax = 2;

T0 = 0 ; T   = 0.5 ;
%%
nstep_x = length( N ) ;
Error_Table = zeros( nstep_x ,5 );
Bound_Table =  zeros( nstep_x ,3 );

for i = 1: nstep_x
    
    nx = N(i);
    ny = N(i);
  
dx = ( xmax - xmin )/ nx ;
dy = ( ymax - ymin ) / ny ;

X = xmin : dx : xmin + nx*dx ;
Y = ymin : dy : ymin + ny*dy ;
x = X( 2: end-1 ) ; % periodic
y = Y( 2: end-1 ) ;
 [ X1, Y1 ] = meshgrid( Y , X ) ; % (nx+1) * (nx+1)
[ X2, Y2 ] = meshgrid( y, x ) ; %  nx *  nx

% dt = dx/( 6*df ) ;
 nt = N_t( i );
% dt = dx*dx;
dt = ( ( T - T0 )/ nt ) ;
 T = T0 + nt*dt ;


%-- Initial ------------
U0 = Exact_sol( X1, Y1, cx, cy, T0, gamma ) ;
U  =  Exact_sol( X1, Y1, cx, cy, T, gamma ) ;
% u_n = U0  ;

%% common part
miu_x=dt/(dx*dx);
miu_y=dt/(dy*dy);
lam_x=dt/dx;
lam_y=dt/dy;
M=1;
% [Mass,u, kxi_plus,yinta, Gplus, Count] = Solve_u_inter_DYak_S( miu_x, miu_y,  lam_x, lam_y, nx , ny,dx, dy,  gamma, nt, dt, X2, Y2, x,y,U0, T0,T) ;
 [Mass_u,u, kxi_plus,yinta, Gplus, Mass_err_u] = new_Solve_u_inter_DYak_S( miu_x, miu_y,  lam_x, lam_y, nx , ny,dx, dy,  gamma, nt, dt, X, Y,x,y, U0, T0,T,M);

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

M_err = max(max(U0))-max(max(u)) ;  %两个结果均为正，则说明 用了限制器 结果会保界
m_err = min(min(u))-min(min(U0)) ;


    % put n in the 1st column of table
    Error_Table(i,  1) = dt;
    % compute the Lmax norm of error  and put it in the 4th column
    Error_Table(i, 2)  =  Emax;
    % compute the L2 norm of error and put it in the 2nd column
    Error_Table(i, 4)  =  E2 ;

    Bound_Table(i,  1) = dt;
    Bound_Table(i, 2) =  M_err;
    Bound_Table( i, 3 ) =  m_err;

    if i ~=1
        Error_Table(i, 3) = log(Error_Table(i, 2)/Error_Table(i-1, 2))...
            /log(Error_Table(i, 1)/Error_Table(i-1, 1));

        % compute the ratio for E2, and put it in the 5th column of table
        Error_Table(i, 5) = log(Error_Table(i, 4)/Error_Table(i-1, 4))...
            /log(Error_Table(i,1)/Error_Table(i-1, 1));

    end


end

format short e
% convert array to table for visualization
TableResults_Err = array2table(Error_Table, 'VariableNames', ...
    {'dt', 'E_2', 'ratio E_2','E_max','ratio E_max'})
%--------------------------------------------------
TableResults_Bound = array2table(Bound_Table, 'VariableNames', {'dt', 'M_err', 'm_err'})





