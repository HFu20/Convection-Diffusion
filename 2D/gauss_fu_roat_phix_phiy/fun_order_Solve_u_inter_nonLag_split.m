function [E2, Emax, wtime, m_err ,M_err] = fun_order_Solve_u_inter_nonLag_split(nx, gamma )

%  nx = 200;
% gamma = 5*1e-4;  
sigma2 =  1e-3 ;  % 2*sigma2 = 1e-3
ny = nx ;  
xmin =-1; xmax =  1 ;
ymin = xmin; ymax = xmax  ;
dx = ( ( xmax - xmin )/ nx )  ; dy =  ( ( ymax - ymin )/ ny )  ;
T0 = 0 ; T   = 2;
   dt =  dx^2  ;
%   dt =  dx/(6*df) ; % dx*5/(2*0.8596)
tt = T0 : dt : T ; nt  = length(tt)-1;
% nt =300;
% dt = (T-T0)/nt ;   tt = T0 : dt : T ;
X = xmin : dx : xmin + nx*dx ;
Y = ymin : dy : ymin + ny*dy ;
x = X( 2: end-1 ) ;  y = Y( 2: end-1 ) ;
% [Y1,X1] = meshgrid( x,y) ;
% [X2,Y2] =  meshgrid( x,y) ;
% x = X( 2: end-1 ) ;  y = Y( 2: end-1 ) ;
 
U0 = Exact_sol( x, y,  T0, gamma ,sigma2 ) ;
U  =  Exact_sol( x, y,  T, gamma,sigma2  ) ; %不是实际解当Source=0时
% U0 = Exact_sol_2(  x, y,cx, cy, T0, gamma ) ;
% U  =  Exact_sol_2( x, y, cx, cy, T, gamma ) ; %不是实际解当Source=0时
u_n  = U0 ;


miu_x=dt/(dx*dx); miu_y=dt/(dy*dy);
lam_x = dt/dx; lam_y = dt/dy ;


tic
%----matrix -------------
% [ matrix_Delta_x_2,matrix_Delta_x ,matrix_A1x ,matrix_A2x ] =  Matrix_periodic(nx) ;
% [matrix_Delta_y_2,matrix_Delta_y ,matrix_A1y ,matrix_A2y ] =  Matrix_periodic(ny) ;
[ matrix_Delta_x_2,matrix_Delta_x ,matrix_A1x ,matrix_A2x ] =  Matrix_interior(nx-1) ;
[matrix_Delta_y_2,matrix_Delta_y ,matrix_A1y ,matrix_A2y ] =  Matrix_interior(ny-1) ;


H_1x = matrix_A1x - gamma*miu_x/4*matrix_Delta_x_2;
H_1y = matrix_A1y - gamma*miu_y/4*matrix_Delta_y_2;
H_2x = matrix_A1x + 1/4*gamma*miu_x*(matrix_Delta_x_2);
 H_2y = matrix_A1y + 1/4*gamma*miu_y*(matrix_Delta_y_2);

 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
% kxi_plus = -1e-12 ;
 
% tt = T0 : dt : T ;
% Mass_u = zeros( length( tt ),1 ); % 初始质量也要计算
%  Mass_S = zeros( length( tt ),1 ); % 初始质量也要计算
%  %初始质量
%  Mass_u(1) =  sum( sum(u_n ) )  *dx*dy ; % 只考虑内点质量（边界为0）
% 
%  for k = 0 : nt
%  Mass_S(1:end) =    sum(sum( fun_Source( X, Y, gamma , T0+k*dt  ,sigma2 )  ) )*dx*dy ;
%  end
% % Mass_S(1) = 0;
% %  yinta只跟t有关，故第k层时，只是一个数
% Mass_S_half = 0.5*dt*( Mass_S(2:end) + Mass_S(1:end-1) ) ;
% %%%%%%%%%%%%%%%%%%%

 
tic
for k = 1:nt

    %保留A_1y*unPlusHf1,第2步里直接运用
%     u_bar_nPlusHf1 = H_2x*u_n*matrix_A1y'+ miu_y *gamma/ 2*matrix_A1x*u_n*matrix_Delta_y_2' ;
%       u_bar_nPlusHf1 =  H_1x\u_bar_nPlusHf1 ;   %没啥效果
%     u_nPlusHf2 = (u_bar_nPlusHf1-miu_y*gamma/4* u_n*matrix_Delta_y_2') / H_1y' ;
    u_bar_nPlusHf1 =  H_1x\( H_2x*u_n*H_2y'  )  ; 
    u_nPlusHf2       = u_bar_nPlusHf1/ H_1y' ;



    %% ---------------------------------------------------------------
    u_nPlusHf3_bar   =  matrix_A2x*u_nPlusHf2*matrix_A2y'-matrix_Delta_x*lam_x* fun_f(u_nPlusHf2,x,y)*matrix_A2y'...
        - lam_y* matrix_A2x*(fun_g(u_nPlusHf2,x,y)) *matrix_Delta_y'  ...
        + dt*matrix_A2x*fun_Source (x,y, gamma, T0+(k-1)*dt ,sigma2  )*matrix_A2y' ;  %g
    u_n3  =  matrix_A2x\ u_nPlusHf3_bar;  %会变快
    u_nPlusHf3   =     u_n3/matrix_A2y' ;

    u_nPlusHf4_bar   = matrix_A2x* u_nPlusHf2*matrix_A2y' - matrix_Delta_x*lam_x* 0.5*(fun_f(u_nPlusHf2,x,y)+fun_f(u_nPlusHf3,x,y) )*matrix_A2y'  ...%f
        -   lam_y/2*  matrix_A2x*(fun_g(u_nPlusHf3,x,y)+fun_g(u_nPlusHf2,x,y) ) * matrix_Delta_y'  ...
        + dt/2 * matrix_A2x*( fun_Source (x,y, gamma, T0+ (k-1)*dt,sigma2   ) + fun_Source (x,y,  gamma, T0+(k)*dt,sigma2   ) )*matrix_A2y';%g;
    u_n4   =    matrix_A2x\ u_nPlusHf4_bar ;
    u_nPlusHf4   =     u_n4 /matrix_A2y' ;
    %% ---------------------------------------------------------------
%     u_bar_nPlusHf5  =   H_2x*u_nPlusHf4 *matrix_A1y'+miu_y*gamma/ 2*matrix_A1x*u_nPlusHf4 *matrix_Delta_y_2';
%    u_bar_nPlusHf5 = H_1x\u_bar_nPlusHf5 ;
%     u_n =  (u_bar_nPlusHf5-miu_y/4*gamma*u_nPlusHf4*matrix_Delta_y_2')/ H_1y' ;

   u_bar_nPlusHf5 =  H_1x\( H_2x*u_nPlusHf4*H_2y'  )    ;
    u_n     = u_bar_nPlusHf5/ H_1y'  ; 
 
%        Mass_u(k+1) =  sum( sum(u_n ) )  *dx*dy ; % 只考虑内点质量（边界为0）
   
end
 
%  Mass_err= Mass_u(2:nt+1) - Mass_u(1:nt )- Mass_S_half ;



wtime = toc;
fprintf (1,'\n Split_nonADI took %f seconds with Nx = %d !!!\n',wtime,nx);
% 
% figure
% plot(tt(2:end),Mass_err,LineWidth=1.2)
% legend('Mass_{err}')
% 
% figure
% surf( x,y,u_n' ) ;
% xlabel("x-axis"); ylabel("y-axis")
% %title("Numerical solution")
% shading interp
% view(2)
%  colorbar
% % ticks = linspace(0, 1, 6); % 5个刻度
% %  colorbar('Ticks', ticks);

figure
surf( x,y,u_n' ) ;
xlabel("x-axis"); ylabel("y-axis")
% title("Numerical solution")
shading interp
colorbar
% 
% figure
% surf( x,y,U0' ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("Initial solution")
% shading interp
% colorbar
 
figure
contour( x,y,u_n' ) ;
xlabel("x-axis"); ylabel("y-axis")
title("Numerical solution")
shading interp
colorbar
% % 
% figure
% contour( x,y, U' ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("Exact solution")
% shading interp
% colorbar

% figure
% surf( x,y,u_n'-U' ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("u-U")
% shading interp
% 
% xx = x(1:nx/20:end) ;
% yy = y(1:ny/20:end) ;
% % ff=  -pi*cos((pi*xx')/2)*sin((pi*yy)/2) ;
% % gg =   pi*sin((pi*xx')/2)*cos((pi*yy)/2) ;
%   [X1,Y1] = meshgrid( xx,yy ) ;
% f =-pi*cos((pi*X1)/2).*sin((pi*Y1)/2)  ;
% g =  pi*sin((pi*X1)/2).*cos((pi*Y1)/2) ;
% figure
% quiver(xx, yy,  f,g);
% % axis equal
% xlabel("x-axis"); ylabel("y-axis")

% figure
% quiver(xx, yy,gg, ff );
% axis equal
%-- Error ------
Emax = max(max(abs(U-u_n))) ;
E2 = sqrt(sum(dy*(sum(dx*((U-u_n).^2))))) ;
m_err = min(min(u_n))-min(min(U0)) ;
M_err = max(max(U0))-max(max(u_n)) ;

fprintf(1,'\n [Emax, E2] = [ %e, %e], m_err = %e,M_err = %e \n',Emax ,E2,m_err ,M_err );

end
