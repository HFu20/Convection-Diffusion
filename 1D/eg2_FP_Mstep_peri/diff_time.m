%--初值为 sin(2*pi*x) 在不同时刻画的图
clear
%------------------
gamma = 1; 
xmin = -2*pi;
xmax = 2*pi ;
cx = 1;

T0 = 0 ;
T=2;

df = 1;
nx= 80;

dx= ((xmax  - xmin)/nx);
%dt=dx/(df);
 nt =40 ;

X=( xmin:dx:xmax)';
x=X(1:end-1); % 
%nx_1=nx-1;

max_u = zeros( 4,1 );
min_u = zeros( 4,1 ) ;
% 
% % [u_1,~] = fun_eg1_figure_gauss(T,dx,dt,xmin,xmax);
% % [u_2,~]= fun_eg1_figure_gauss(T/2,dx,dt,xmin,xmax);
% % [u_3,U0] = fun_eg1_figure_gauss(T/4,dx,dt,xmin,xmax);
% [u_1 ,U0, Mass_1,tt_1,  max_u(1), min_u(1)]       = fun_perio_FP_U0_sin1( x ,  nx, dx, nt,0.01,T0, gamma ); % Generate n*n matrix
% [u_2 ,~,  Mass_2,  tt_2,  max_u(2), min_u(2)]  =  fun_perio_FP_U0_sin1(  x ,  nx, dx, nt,0.1,T0, gamma );
% [u_3 ,~,  Mass_3, tt_3 ,  max_u(3), min_u(3)]   =   fun_perio_FP_U0_sin1(  x ,  nx, dx, nt,0.4,T0, gamma );
% % [u_4 ,~, Mass_4 , tt_4, max_u(4), min_u(4)] =  fun_perio_U0_sin1( x ,  nx_1, dx, dt,T*4, T0,gamma );

% figure
% % plot(X,U0,'k')
% % hold on
% %plot(x,U,'b');
%  hold on
% plot(x,u_1,'k', 'linewidth', 1 );
% plot(x,u_2,'r-d', 'linewidth', 1 );
% plot(x,u_3,'b-*', 'linewidth', 1 );
% % plot(x,u_4,'-')
% xlim([xmin xmax])
% xlabel('x')
% ylabel('u')
% legend('T=0.01',...
%     ' T=0.1',...
%     ' T=0.4')
% 
%  figure
% plot(tt_3, Mass_3, '-')
% ylabel('Mass_{err}') %已达到机器误差
%  ylim([Mass_3(1)-1e-12 Mass_3(1)+1e-12])
% %yticks( [-1e-12,-0.5*1e-12,0,0.5*1e-12,1e-12 ] )
%  xlabel('Time')

%% 比较加乘子前后的质量
[uu_1 ,new_Mass_1,tt,  max_u(1), min_u(1), count_1,Energy_1,dt]      = fun_perio_FP_Lag( x ,  nx, dx, nt,T,T0, gamma ); 
[uu_2, new_Mass_2, ~,  max_u(2), min_u(2), count_2, Energy_2]       = fun_perio_FP_non( x ,  nx, dx, nt,T,T0, gamma ); 

figure
plot(x, uu_1,'r-d','Markersize',8)
hold on
plot(x, uu_2,'b-o','Markersize',5)
xlim([-2*pi 2*pi])
xlabel('x-axis')
%ylabel('u')
legend('With multipliers','Without multipliers')

 figure
plot(tt, new_Mass_1-new_Mass_1(1), 'g-^',LineWidth=1.2)
hold on 
plot(tt, new_Mass_2-new_Mass_2(1), 'b-o',LineWidth=1) 
% ylim([-5e-13  5e-13])
 xlabel('Time')
ylabel('Mass_{err}')
legend('With multipliers','Without multipliers')

figure
% E_err =  Energy   - Energy(1) ;
plot( tt,   real(Energy_1),'r-.',  LineWidth=1.2)
hold on 
plot( tt,   real(Energy_2),'b--' ,LineWidth=1.2)
ylabel('Energy') %已达到机器误差
%  ylim([-1e-14  1e-14])
% yticks( [-1e-10,-0.5*1e-10,0,0.5*1e-10,1e-10 ] )
xlim([ dt,  T])
xlabel('Time')
legend('With multipliers','Without multipliers')