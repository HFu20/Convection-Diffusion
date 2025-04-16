% 与沈捷老师的BDF2比较

 function [u ] = fun_bdf2_2d( miu_x, miu_y,  lam_x, lam_y, nx , ny, gamma, nt, dt, X2, Y2,U0, T0) 
%%

% -- Solve ----------
tic

 %----matrix -------------
[ matrix_Delta_x_2,matrix_Delta_x ,matrix_A12x ,matrix_A6x ] = Matrix_interior( nx-1 );
[ matrix_Delta_y_2,matrix_Delta_y ,matrix_A12y ,matrix_A6y ] = Matrix_interior( ny-1 );
 u  = U0( 2:end-1, 2:end-1 ) ;
 uminus = u(:) ;
 Iy = eye( ny-1,ny-1 );
 Ix = eye( nx-1, nx-1 ) ;

 %% BDF1 %%%%%%%%%%%%%%%%%
A_all = kron( Iy, matrix_A6x*matrix_A12x )*kron(  matrix_A6y*matrix_A12y, Ix ) ;
 Left_bdf1 = A_all  - miu_x*gamma*kron( matrix_A6y*matrix_A12y, Ix ) *kron( Iy, matrix_A6x*matrix_Delta_x_2 )...
                       - miu_y*gamma*kron( Iy, matrix_A6x*matrix_A12x )*kron( matrix_A6y*matrix_Delta_y_2,Ix ) ;
i = 1 ;
 S =  fun_Source( X2, Y2, gamma ,T0+ i*dt ) ;
 un = Left_bdf1\(  A_all*uminus ...
                       -lam_x*kron(  matrix_A6y*matrix_A12y, Ix )*kron( Iy, matrix_A12x*matrix_Delta_x )* fun_f(uminus)...
                       -lam_y*kron( Iy, matrix_A6x*matrix_A12x )* kron( matrix_A12y*matrix_Delta_y, Ix )*fun_g(uminus) ...
                       + dt* A_all*S(:) );
 
 Left_bdf2 = 3/2*A_all - miu_x*gamma*kron( matrix_A6y*matrix_A12y, Ix ) *kron( Iy, matrix_A6x*matrix_Delta_x_2 )...
                        - miu_y*gamma*kron( Iy, matrix_A6x*matrix_A12x )*kron( matrix_A6y*matrix_Delta_y_2,Ix ) ;
 for i = 2 : nt
  S =  fun_Source( X2, Y2, gamma ,T0+ i*dt ) ;
  ubar = 2*un - uminus ;
   uplus = Left_bdf2\( 2*A_all*un - 1/2*A_all*uminus ...
                 -  lam_x*kron(  matrix_A6y*matrix_A12y, Ix )*kron( Iy, matrix_A12x*matrix_Delta_x )* fun_f(ubar) ...
                 - lam_y*kron( Iy, matrix_A6x*matrix_A12x )* kron( matrix_A12y*matrix_Delta_y, Ix )*fun_g(ubar) ...
                 + dt*A_all*S(:) ) ;
   uminus = un ;
   un = uplus;
 end
toc
uu = reshape( uplus, [nx-1, ny-1] ) ;
u = U0;
u(2:end-1, 2:end-1) = uu ;
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


%-- Error ------
% Emax = max(max(abs(U-u)));
% E2 = sqrt(sum(dy*(sum(dx*((U-u).^2))))) ;
% M_err = max(max(U0))-max(max(u))  %两个结果均为正，则说明 用了限制器 结果会保界
% m_err = min(min(u))-min(min(U0))
