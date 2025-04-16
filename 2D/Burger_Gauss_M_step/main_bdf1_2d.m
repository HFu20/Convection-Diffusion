% 与沈捷老师的BDF2比较
clear

%%
% gamma =0.005;
gamma =1;
df = 1 ;
cx = 1 ;
cy = 1 ;

xmin = -1 ;
xmax = 2 ;
ymin = -1 ;
ymax = 2;

T0 = 0 ;
T   = 0.1 ;

nx = 300 ;
ny = 300 ;
dx = ( xmax - xmin )/ nx ;
dy = ( ymax - ymin ) / ny ;

%按推导的步长比
% dx = gamma/(5*df) ;
% dy =gamma/(5*df) ;
% nx =  floor( ( xmax - xmin )/ dx );
% ny =  floor( ( ymax - ymin )/ dy );


X = xmin : dx : xmin + nx*dx ;
Y = ymin : dy : ymin + ny*dy ;
x = X(2 : end-1 ) ;  
y = Y( 2: end-1 ) ;
[ X1, Y1 ] = meshgrid( Y , X ) ; % (nx+1) * (nx+1)
[ X2, Y2 ] = meshgrid( y, x ) ; %  nx *  nx

% dt = min( dx/( 3*df ) , 5/(6*gamma) *dx^2  ) ;
%  dt = dx/(6* df ) ;
   dt = dx^2 ;
nt = floor( ( T - T0 )/ dt ) ;
T = T0 + nt*dt ;
tt = T0 : dt : T ;

%-- Initial ------------
U0 = Exact_sol( X1, Y1, cx, cy, T0, gamma ) ; 
U  =  Exact_sol( X1, Y1, cx, cy, T, gamma ) ; %不是实际解当Source=0时
%u_n = U0  ;

Mass_U = zeros(nt+1,1) ;
Mass_S = zeros(nt+1,1) ;
for k = 1 : nt+1
    t =  T0+ (k-1)*dt ;
    U  =  Exact_sol( X1, Y1, cx, cy,t, gamma ) ; 
for i =1 : length(U(:,1))
    for j = 1 : length(U(1,:))
        Mass_U(k) = Mass_U(k) + U(i,j)*dx*dy ; % 只考虑内点质量（边界为0）
       %  Mass_S(k) =  Mass_S(k) + fun_Source(X(i), Y(j), gamma , T0 )*dt*dx*dy ;       
    end
end
end
Mass_exc = Mass_U - Mass_S ;

% figure
% plot(  tt, Mass_U-Mass_U(1) );
% legend('Exact Mass_{err}') %已达到机器误差
%  ylim([ -1e-12 1e-12])
% %   yticks( [ Mass(1)-1e-12 ,0,Mass(1)+1e-12 ] )
%  xticks( T0: (T-T0)/4 : T  )
%  xlabel('Time')

%% common part
miu_x=1/(dx*dx);
miu_y=1/(dy*dy);
lam_x=1/dx;
lam_y=1/dy;


% -- Solve ----------
tic
 % BDF1
 %----matrix -------------
[ matrix_Delta_x_2,matrix_Delta_x ,matrix_A12x ,matrix_A6x ] = Matrix_interior( nx-1 );
[ matrix_Delta_y_2,matrix_Delta_y ,matrix_A12y ,matrix_A6y ] = Matrix_interior( ny-1 );
 u  = U0( 2:end-1, 2:end-1 ) ;
 u = u(:) ;
 Iy = eye( ny-1,ny-1 );
 Ix = eye( nx-1, nx-1 ) ;
 Left_1 = kron( Iy, matrix_A6x*matrix_A12x )*kron(  matrix_A6y*matrix_A12y, Ix )  ...
                - dt*gamma*miu_x*kron( matrix_A6y*matrix_A12y, Ix ) *kron( Iy, matrix_A6x*matrix_Delta_x_2 )...
                - dt*gamma*miu_y*kron( Iy, matrix_A6x*matrix_A12x )*kron( matrix_A6y*matrix_Delta_y_2,Ix ) ;

 for i =1 : nt
   f_u =  fun_f(u) ;
   g_u = fun_g(u) ;
   S =  fun_Source( X2, Y2, gamma ,T0+ i*dt ) ;
   u = Left_1\(  kron( Iy, matrix_A6x*matrix_A12x )*kron(  matrix_A6y*matrix_A12y, Ix )*u ...
                       -dt*lam_x*kron(  matrix_A6y*matrix_A12y, Ix )*kron( Iy, matrix_A12x*matrix_Delta_x )*f_u...
                       -dt*lam_y*kron( Iy, matrix_A6x*matrix_A12x )* kron( matrix_A12y*matrix_Delta_y, Ix )*g_u ...
                       + dt* kron( Iy, matrix_A6x*matrix_A12x )*kron(  matrix_A6y*matrix_A12y, Ix )*S(:) );
 end
toc
 un = reshape( u, [nx-1, ny-1] ) ;
u = U0;
u(2:end-1, 2:end-1) = un ;
% -- Figure -----------------
% note X Y

figure
surf( X1, Y1, U0 ) ;
xlabel("x-axis"); ylabel("y-axis")
title("Initial solution")
shading interp

figure
surf( X1, Y1, u ) ;
xlabel("x-axis"); ylabel("y-axis")
title("Numerical solution")
shading interp

figure
surf( X1, Y1, U ) ;
xlabel("x-axis"); ylabel("y-axis")
title("Exact solution")
shading interp

figure
surf( X1, Y1, U-u ) ;
xlabel("x-axis"); ylabel("y-axis")
title("Exact solution")
shading interp

%  figure
%  surf( X1, Y1, min(u,0) ) ;
%  xlabel("x-axis"); ylabel("y-axis")
%  title("min(u,0)")
%  colorbar
%  shading interp

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

%-- Error ------
Emax = max(max(abs(U-u)))
E2 = sqrt(sum(dy*(sum(dx*((U-u).^2)))))
E3 = max(  max( abs(U-u)/max( abs( U )) )  ) 
% M_err = max(max(U0))-max(max(u))  %两个结果均为正，则说明 用了限制器 结果会保界
% m_err = min(min(u))-min(min(U0))
