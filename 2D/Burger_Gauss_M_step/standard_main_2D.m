
clear

%%
%23.04.03
%limter1中 T1(j)<lbound-1e-20 eps在小就不符合
%  加限制器仍不保界 甚至Non
% 23.04.05
%取源项为0 
% 在eps条件下 舍掉低于机器精度的值
%%
% gamma =0.005;
gamma =0.005;
df = 0.4496 ;
cx = 1 ;
cy = 1 ;
niu = 5/( 6*gamma ) ;

xmin = -1 ;
xmax = 2  ;
ymin = -1 ;
ymax = 2  ;

T0 = 0 ;
T   = 0.6 ;

nx = 200 ;
 ny = nx ;
dx = ( xmax - xmin )/ nx ;
dy = ( ymax - ymin ) / ny ;

% dt = min( dx/( 3*df ) , 5/(6*gamma) *dx^2  ) ;
M = 1  ;
dt =  dx/(6*df ) ;
%    dt = dx^2 ;
nt = floor( ( T - T0 )/ dt ) ;
T = T0 + nt*dt ;
tt = T0 : dt : T ;

%按推导的步长比
% dx = gamma/(5*df) ;
% dy =gamma/(5*df) ;
% nx =  floor( ( xmax - xmin )/ dx );
% ny =  floor( ( ymax - ymin )/ dy );


X = xmin : dx : xmin + nx*dx ;
Y = ymin : dy : ymin + ny*dy ;
x = X(1 : end-1 ) ; % periodic
y = Y( 1: end-1 ) ;
%[ X1, Y1 ] = meshgrid( Y , X ) ; % (nx+1) * (nx+1)
[ X2, Y2 ] = meshgrid( y, x ) ; %  nx *  nx



%-- Initial ------------
U0 = Exact_sol( X2, Y2, cx, cy, T0, gamma ) ; 
U  =  Exact_sol( X2, Y2, cx, cy, T, gamma ) ; %不是实际解当Source=0时
 
% figure
% surf( X2, Y2, U ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("Exact solution")
% shading interp
% 
Mass_U = zeros(nt+1,1) ;
Mass_S_ex = zeros(nt+1,1) ;
for k = 1 : nt+1
    t =  T0+ (k-1)*dt ;
    U  =  Exact_sol( X2, Y2, cx, cy,t, gamma ) ; 
%  % Mass_U_Left  = (U^{n+1},1 );
%  % Mass_U_Right  = (U^{n},1 ) + dt*1/2*(fun_S^{n+1}+fun_S^n );
%  %理论上，应有 Mass_U_Left = Mass_U_Right
    Mass_U(k) = sum( sum(U ) )*dx*dy ; % 只考虑内点质量（边界为0）
     Mass_S_ex(k )    =   sum( sum(fun_Source( X2, Y2, gamma , t  )  ) ) *dx*dy ;   
end
% %  
  Mass_err_U = Mass_U(2:nt+1) - Mass_U(1:nt )- dt* ( Mass_S_ex(2:nt+1) + Mass_S_ex(1:nt))/2 ;
% % 
figure
plot(  tt(2:end), Mass_err_U );
legend('Exact Mass_{err}') %已达到机器误差
%  ylim([ -1e-16 1e-16])
%   yticks( [ Mass(1)-1e-12 ,0,Mass(1)+1e-12 ] )
 xticks( T0: (T-T0)/4 : T  )
 xlabel('Time')

 figure
plot(  tt, Mass_U );
legend('Exact Mass_{U}') %已达到机器误差
%  ylim([ -1e-14 1e-14])
%   yticks( [ Mass(1)-1e-12 ,0,Mass(1)+1e-12 ] )
 xticks( T0: (T-T0)/4 : T  )
 xlabel('Time')

%% common part
miu_x=dt/(dx*dx);
miu_y=dt/(dy*dy);
lam_x=dt/dx;
lam_y=dt/dy;


% -- Solve ----------
tic
%   [Mass,u, kxi_plus,yinta, Gplus, Count] = Solve_u_perio_DYak_S( miu_x, miu_y,  lam_x, lam_y, nx , ny,dx, dy,  gamma, nt, dt, X2,Y2, U0, T0,T);
  [Mass_u,u, kxi_plus,yinta, Gplus,Mass_err_u] = mass_Solve_u_inter_DYak_S( miu_x, miu_y,  lam_x, lam_y, nx , ny,dx, dy,  gamma, nt, dt, X2, Y2, x,y,U0, T0,T,M) ;
%  [u ] = u_Solve_u_inter_DYak_S( miu_x, miu_y,  lam_x, lam_y, nx , ny,   gamma, nt, dt, X2, Y2,  U0, T0, M) ;
toc

% for i = 1 : nx 
%     for j =1:ny  
%      if u(i,j)> max(max(U0))
%          u(i,j)= max(max(U0)) ;
%      elseif u(i,j)< min(min(U0))
%          u(i,j) = min(min(U0)) ;
%      end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%
% -- Figure -----------------
% note X Y

% figure
% surf( X1, Y1, U0 ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("Initial solution")
% shading interp

figure
surf( X2, Y2, u ) ;
xlabel("x-axis"); ylabel("y-axis")
%title("Numerical solution")
shading interp

% figure
% contour( X2, Y2, u ) ;
% xlabel("x-axis"); ylabel("y-axis")
% %title("Numerical solution")
% shading interp
% colorbar
% 
% 
figure
plot(  tt(2:end), Mass_err_u)
xlabel("Time")
legend('Numerical Mass_{err}') %已达到机器误差
 % ylim([ -5e-14 5e-14])

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


%-- Error ------

fprintf('nx = %d, T = %f \n',nx,T )
Emax = max(max(abs(U-u))) ;
E2 = sqrt(sum(dy*(sum(dx*((U-u).^2))))) ;
% M_err = max(max(U0))-max(max(u))  %两个结果均为正，则说明 用了限制器 结果会保界
m_err = min(min(u))-min(min(U0)) ;
fprintf(1,'\n [Emax, E2] = [ %lo, %lo]\n',Emax ,E2  );

