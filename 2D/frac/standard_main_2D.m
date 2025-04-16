clear
%---data -----------------
gamma =0.01;
df = 1 ;
cx = 1 ;
cy = 1 ;
niu = 5/( 6*gamma ) ;

xmin = -1.5 ;
xmax = 1.5 ;
ymin = -1.5 ;
ymax = 1.5;
% xmin = -3 ;
% xmax = 3 ;
% ymin = -3 ;
% ymax = 3;

T0 = 0 ;
T   = 0.5;

nx =400 ;
ny = nx ;
dx = ( xmax - xmin )/ nx ;
dy = ( ymax - ymin ) / ny ;
% dx = 1/128 ;
% dy = 1/128 ;
% nx = floor( ( xmax - xmin )/ dx ) ;
% ny = floor( ( ymax - ymin ) / dy  );
%按推导的步长比
% dx = gamma/(5*df) ;
% dy =gamma/(5*df) ;
% nx =  floor( ( xmax - xmin )/ dx );
% ny =  floor( ( ymax - ymin )/ dy );


X = xmin : dx : xmin + nx*dx ;
Y = ymin : dy : ymin + ny*dy ;
x = X(2 : end-1 ) ; % periodic
y = Y( 2: end-1 ) ;
[ X1, Y1 ] = meshgrid( Y , X ) ; % (nx+1) * (nx+1)
[ X2, Y2 ] = meshgrid( y, x ) ; %  nx *  nx

% dt = min( dx/( 3*df ) , 5/(6*gamma) *dx^2  ) ;
K = 1;
dt = dx/(3*df ) ;
% dt = 0.001 ;
nt = floor( ( T - T0 )/ dt ) ;
T = T0 + nt*dt ;

%-- Initial ------------
U0 = zeros( nx+1, ny+1 ) ;
for i =1 : nx+1 
    for j =1: ny+1
if X(i)^2 + Y(j)^2 < 0.5
    U0(i,j) = 1;
else
    U0(i,j) = 0;
end
    end
end

%% common part
miu_x=dt/(dx*dx);
miu_y=dt/(dy*dy);
lam_x=dt/dx;
lam_y=dt/dy;

tt = T0 : dt : T ;
% -- Solve ----------
tic
   [Mass,u, kxi_plus,yinta, Gplus, Count] = Solve_u_inter_DYak_S( miu_x, miu_y,  lam_x, lam_y, nx , ny,dx, dy,  gamma, nt, dt, X2, Y2, U0, T0,T,K) ;
 %                [Mass,u ] = Solve_u_inter_non( miu_x, miu_y,  lam_x, lam_y, nx , ny,dx, dy,  gamma, nt, dt, X2, Y2, U0, T0,T);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%
% eps = 1e-10;
% for i =1 : nx
%     for j = 1 : ny
%         if ( u(i,j) <= max( max(U0) )  + eps)  && ( u(i,j) > max( max(U0) )  )
%             u(i,j) = max( max(U0) ) ;
%         elseif  ( u(i,j) >=  min( min(U0) ) -eps ) && ( u(i,j) < min( min(U0) )  )
%             u(i,j) = min( min(U0) ) ;
%         end
% 
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%
% -- Figure -----------------
% note X Y
% 
% figure
% surf( X1, Y1, U0 ) ;
% xlabel("x-axis"); ylabel("y-axis")
% ylim([ymin ymax])
% xlim([xmin xmax])
% shading interp
%  colorbar
% 
% figure
% contour( X1, Y1, U0 ) ;
% xlabel("x-axis"); ylabel("y-axis")
% ylim([ymin ymax])
% xlim([xmin xmax])
% shading interp
%  colorbar

figure
contour( X1, Y1, u, 20)
 colorbar
%contour( X1, Y1, u ) ;
xlabel("x-axis"); ylabel("y-axis")
%title("Numerical solution")
%shading interp
ylim([ymin ymax])
xlim([xmin xmax])
%saveas(gcf,'u')
 colorbar

figure
mesh( X1, Y1, u ) ;
xlabel("x-axis"); ylabel("y-axis")
%title("Numerical solution")
 colorbar
%shading interp
ylim([ymin ymax])
xlim([xmin xmax])
%saveas(gcf,'u')

%  figure
%  surf( X1, Y1, min(u,zeros(nx+1,ny+1)) ) ;
%  xlabel("x-axis"); ylabel("y-axis")
%  title("min(u,0)")
%  colorbar
%  shading interp
%saveas(gcf,'min(u,0)')

% figure
% plot(  tt,  Mass )
%  xticks( T0: (T-T0)/4 : T  )
%  xlabel('Time') 
% title('Mass') %已达到机器误差
%  ylim([Mass(1)-1e-8  Mass(1)+1e-8])
% %   yticks( [ Mass(1)-1e-12 ,0,Mass(1)+1e-12 ] )
%  saveas(gcf,'Mass')
% % 
% figure
% plot(  tt,  Mass-Mass(1) );
% ylabel('Mass_{err}') %已达到机器误差
% xlabel('Time')
% xticks( T0: (T-T0)/4 : T  )
% ylim([-1e-8  1e-8])
%  yticks( [ -1e-8 ,0, 1e-8 ] )



%-- Error ------
% Emax = max(max(abs(U-u)))
% E2 = sqrt(sum(dy*(sum(dx*((U-u).^2)))))
% E3 = max(  max( abs(U-u)/max( abs( U )) )  ) 
M_err = max(max(U0))-max(max(u))  %两个结果均为正，则说明 用了限制器 结果会保界
m_err = min(min(u))-min(min(U0))
