% 与沈捷老师的BDF2比较
 clear
% function []=mu
%%
% gamma =0.005;
gamma =0.005;
df =  0.4496 ;
cx = 1 ;
cy = 1 ;
niu = 5/( 6*gamma ) ;

xmin = -1 ;
xmax = 2 ;
ymin = -1 ;
ymax = 2;

T0 = 0 ;
T   = 0.6 ;

nx = 200 ;
ny = nx ;
dx = ( xmax - xmin )/ nx ;
dy = ( ymax - ymin ) / ny ;

X = xmin : dx : xmin + nx*dx ;
Y = ymin : dy : ymin + ny*dy ;
x = X(1 : end-1 ) ;  
y = Y( 1: end-1 ) ;
%[ X1, Y1 ] = meshgrid( Y , X ) ; % (nx+1) * (nx+1)
[ X2, Y2 ] = meshgrid( y, x ) ; %  nx *  nx

% dt = min( dx/( 3*df ) , 5/(6*gamma) *dx^2  ) ;
     dt = dx/(6* df ) ;
%   dt = dx^2 ;
nt = floor( ( T - T0 )/ dt ) ;
T = T0 + nt*dt ;
tt = T0 : dt : T ;

%-- Initial ------------
U0 = Exact_sol( X2, Y2, cx, cy, T0, gamma ) ; 
U  =  Exact_sol( X2, Y2, cx, cy, T, gamma ) ; %不是实际解当Source=0时
 
%% common part
miu_x=dt/(dx*dx);
miu_y=dt/(dy*dy);
lam_x=dt/dx;
lam_y=dt/dy;


% -- Solve ----------
tic

 %----matrix -------------
[ matrix_Delta_x_2,matrix_Delta_x ,matrix_A12x ,matrix_A6x ] = Matrix_periodic(nx) ;
[ matrix_Delta_y_2,matrix_Delta_y ,matrix_A12y ,matrix_A6y ] = Matrix_periodic(ny) ;
 u  = U0 ;
 uminus = u(:) ;
 Iy = eye( ny,ny );
 Ix = eye( nx, nx) ;
% Mass_u(1) = sum(uminus)*dy*dx ;
 %% BDF1 %%%%%%%%%%%%%%%%%
A_all = kron( Iy, matrix_A6x*matrix_A12x )*kron(  matrix_A6y*matrix_A12y, Ix ) ;
 Left_bdf1 = A_all  - miu_x*gamma*kron( matrix_A6y*matrix_A12y, Ix ) *kron( Iy, matrix_A6x*matrix_Delta_x_2 )...
                       - miu_y*gamma*kron( Iy, matrix_A6x*matrix_A12x )*kron( matrix_A6y*matrix_Delta_y_2,Ix ) ;
k = 1 ;
 S =  fun_Source( X2, Y2, gamma ,T0+ k*dt ) ;
 un = Left_bdf1\(  A_all*uminus ...
                       -lam_x*kron(  matrix_A6y*matrix_A12y, Ix )*kron( Iy, matrix_A12x*matrix_Delta_x )* fun_f(uminus)...
                       -lam_y*kron( Iy, matrix_A6x*matrix_A12x )* kron( matrix_A12y*matrix_Delta_y, Ix )*fun_g(uminus) ...
                       + dt* A_all*S(:) );
%   Mass_u(k+1) =  sum(un)*dx*dy ;

 Left_bdf2 = 3/2*A_all - miu_x*gamma*kron( matrix_A6y*matrix_A12y, Ix ) *kron( Iy, matrix_A6x*matrix_Delta_x_2 )...
                        - miu_y*gamma*kron( Iy, matrix_A6x*matrix_A12x )*kron( matrix_A6y*matrix_Delta_y_2,Ix ) ;
 for k= 2 : nt
  S =  fun_Source( X2, Y2, gamma ,T0+ k*dt ) ;
  ubar = 2*un - uminus ;
   uplus = Left_bdf2\( 2*A_all*un - 1/2*A_all*uminus ...
                 -  lam_x*kron(  matrix_A6y*matrix_A12y, Ix )*kron( Iy, matrix_A12x*matrix_Delta_x )* fun_f(ubar) ...
                 - lam_y*kron( Iy, matrix_A6x*matrix_A12x )* kron( matrix_A12y*matrix_Delta_y, Ix )*fun_g(ubar) ...
                 + dt*A_all*S(:) ) ;
   uminus = un ;
   un = uplus;
%     Mass_u(k+1) =  sum(un)*dx*dy ;
 end
toc
u = reshape( uplus, [nx, ny] ) ;

% -- Figure -----------------
% note X Y

% figure
% surf( X1, Y1, U0 ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("Initial solution")
% shading interp
% 
% figure
% surf( X1, Y1, u ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("Numerical solution")
% shading interp
% 
% figure
% surf( X1, Y1, U ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("Exact solution")
% shading interp
% 
% figure
% surf( X1, Y1, U-u ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("U-u")
% shading interp

%  figure
%  surf( X1, Y1, min(u,0) ) ;
%  xlabel("x-axis"); ylabel("y-axis")
%  title("min(u,0)")
%  colorbar
%  shading interp
% 
%  figure
%  surf( X1, Y1, min(U,0) ) ;
%  xlabel("x-axis"); ylabel("y-axis")
%  title("min(U,0)")
%  colorbar
%  shading interp

% figure
% plot(  tt,  Mass );
% title('Mass') %已达到机器误差
% ylim([Mass(1)-1e-8  Mass(1)+1e-8])
%   yticks( [ Mass(1)-1e-12 ,0,Mass(1)+1e-12 ] )
%  xticks( T0: (T-T0)/4 : T  )
%  xlabel('Time')
 
% 
% figure
% plot(  tt,  Mass-Mass(1) );
% ylabel('Mass_{err}') %已达到机器误差
%  ylim([ -1e-12 1e-12])
% % ylim([-1e-8  1e-8])
% %  yticks( [ -1e-8 ,0, 1e-8 ] )
% xticks( T0: (T-T0)/4 : T  )
% xlabel('Time')

wtime = toc;
fprintf (1,'\n BDF took %f seconds to run !!!\n',wtime);

%-- Error ------
Emax = max(max(abs(U-u))) ;
E2 = sqrt(sum(dy*(sum(dx*((U-u).^2))))) ;
fprintf(1,'\n [Emax, E2] = [ %lo, %lo]\n',Emax ,E2  );
% M_err = max(max(U0))-max(max(u))  %两个结果均为正，则说明 用了限制器 结果会保界
% m_err = min(min(u))-min(min(U0))
