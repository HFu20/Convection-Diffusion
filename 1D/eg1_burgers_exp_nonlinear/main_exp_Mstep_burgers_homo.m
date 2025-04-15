%方程本身无质量守恒性
clear
 
%nx = 250; gamma = 0.0005;cx =1 ;
   nx =600;  gamma =5e-4; cx =1 ;
%  nx = 1200; gamma = 2e-4;cx =1 ;

xmin = 0;
xmax = 3;
dx=((xmax  - xmin)/nx);
X=( xmin:dx:xmax)';
x=X(2:end-1); % homogeous boundary

df = 0.25;
%T0=1;T=13;
T0=1;T=4; 
% M =1 ; nt =600 ;
% dt =  (T-T0)/nt ;
%      M =1 ; dt = 4.45*dx/( 5* df);
 %     M =10 ; dt = 3* dx/( df);
 M =1 ; dt =  1.58* dx/( 5*  df) ;
% dt = 1e-5;
%    dt = dx^2;
 
nt=floor((T-T0)/dt);
t=(T0:dt:T0+dt*nt);
T = T0+dt*nt ;

%initial and Exact condition
[U0] = Exact_exp(t(1),  X,  gamma) ;
[U] =  Exact_exp(t(end) , X, gamma) ;
u_n=U0( 2 : end-1 ) ;
kxi = 0;
lam_n = zeros( length(u_n), 1);
tt = T0 : dt : T ;
Mass_u = zeros( length( tt ),1 ); % 初始质量也要计算
Mass_U = zeros( length( tt ),1 );
Iter = zeros( length( tt ) -1 ,1 );  % 一共nt次

a = min(U0) ;
b = max(U0) ;
gprime = @(x) (-2*x +a +b) ;
 g = @(x) (b-x).*(x-a) ;

for k =1 :nt+1
  t =  T0+ (k-1)*dt ;
    U  = Exact_exp(t , X, gamma) ;
 Mass_U(k) = sum( U )*dx; 
end
% figure
%  plot(tt,Mass_U)
figure
 plot(tt,Mass_U-Mass_U(1))


for j  = 1 : length(U0) 
    Mass_u(1) = Mass_u(1) + U0(j)*dx ;
end

[ matrix_Hl , matrix_Hr , matrix_Delta_x_2 , matrix_Delta_x, matrix_A12,matrix_A6] = Matrix_nonp(nx-1, gamma, dx, dt) ;

epsilon = 1 ;
T_lag = tic ;
for k = 1:nt

    u_nPlusHf = matrix_Hl\( matrix_Hr * u_n  ); %+ 0.5*matrix_A12*dt*(kxi + lam_n.*gprime(u_n)) ) ;   % 扩散

% % %  [ u_nPlus1_tilde] =...
%         fun_Mstep_convec( u_nPlusHf, dt, matrix_Delta_x,matrix_A6, kxi,T0, cx,x,gamma ,M,k,lam_n,gprime,u_n) ; %分成m个子步
        [ u_nPlus1_tilde] =  fun_Mstep_convec( u_nPlusHf, dt, matrix_Delta_x,matrix_A6, kxi,T0, cx,x,gamma ,M,k) ;
 %      [u_nPlus1_tilde] = TVB_fun_Mstep_convec( u_nPlusHf, dt, matrix_A6, M, df, dx,nx,epsilon  ) ;
  
   u_nPlus1 =matrix_Hl\(matrix_Hr * u_nPlus1_tilde ); % + 0.5*matrix_A12*dt*(kxi + lam_n.*gprime(u_n)) );  

% [ kxi_plus(k),count(k) ] = lam_fun_kxi_plus( dt,  u_nPlus1,u_n,U0, kxi, dx,  lam_n  ) ;
% yinta =  1/1*dt*(  kxi_plus(k) - kxi -lam_n.*g(u_n) ) ;    
%      z = u_nPlus1 + yinta  ;
% kxi = kxi_plus(k) ;
%     for i = 1: length( u_nPlus1 )
%         if z(i) <= min(U0)
%             u_nPlus1 (i) = min(U0);
%               lam_n(i) = 1*(a - (u_nPlus1(i) + yinta(i)))/(dt*gprime(a))  ;
%         elseif z(i) >= max(U0)
%             u_nPlus1 (i) = max(U0);
%              lam_n(i) =1* (b - (u_nPlus1(i) + yinta(i)))/(dt*gprime(b))  ;
%         else
%             u_nPlus1 (i) = z(i) ;
%                lam_n(i) = 0 ;
%         end
%     end
%  Mass_u(k+1) =sum(u_nPlus1)*dx ;

    %--------------------------------------------


    u_n = u_nPlus1;


end


% 
     for i = 1: length( u_nPlus1 ) 
        if u_nPlus1(i) <= min(U0)
            u_nPlus1 (i) = min(U0);  
         elseif u_nPlus1 (i) >= max(U0)
            u_nPlus1 (i) = max(U0);
        end  
     end
      u_n = u_nPlus1;

u = U ;
u( 2:end-1 ) = u_n;
 toc(T_lag) ;
 min(u);
% subplot(1,4,1)
% plot(x,u,'b');
% title("最终时刻的数值解")
% subplot(1,4,2)
% plot(x,U,'r');
% title("最终时刻的精确解")
% subplot(1,4,3)
% plot(x,U-u,'*-');
% subplot(1,4,4)

figure
plot( X ,U,'r--','LineWidth',1.3);
hold on
plot( X ,u,'b-')
legend('Exact soulution u ','Numerical soulution U')
 xlabel('x-axis')

figure
plot( X ,U-u,'b--','LineWidth',1.3);

%  figure
%  Mass_err_u = Mass_u -Mass_u(1) ;
% plot(tt,Mass_err_u)
%  xlabel('Time')
%%'Position', [0.1 0.1 0.8 0.8]
%第1个参数表示图坐标区域左下角至图形窗口左边的占比
%第2个参数表示图坐标区域左下角至图形窗口下边的占比
%第3个参数表示图坐标区域的宽占比
%第4个参数表示图坐标区域的高占比
% figure
% ax1 = axes('Position', [0.1 0.1 0.8 0.8]);
%  axes(ax1) %大窗口绘图
% plot( X,U,'r-.','lineWidth',2)
% hold on 
% plot(X,u,'b--');
% xlabel(' x-axis')
% legend('Exact solution','Numerical solution')
% %ylabel('$$ f(x) $$', 'interpreter', 'latex')
% %title('局部放大图', 'fontsize', 12)
% % grid on
% % x11 = [0.45, 0.3];
% % y11 = [0.45, 0.3];
% % annotation('arrow',x11, y11);
%------------
% ax2 = axes('Position',[0.5 0.34 0.3 0.3]);
% axes(ax2) %小窗口绘图
% a = 193 : 204;
%  x1 = X(a) ;
%  plot( X,U,'r--','lineWidth',2)
% hold on 
% plot( x1, u(a),'b-o' )
% u1 = min( u(a) ) ;u2 = max( u(a) ) ;
% xlim( [x1(1) x1(end)])
% ylim( [u1, u2])
% 


 
% 
% figure
% plot( X ,U,'b-');
% hold on
% plot( X ,u,'ro')
% legend('Exact soulution u ','Numerical soulution U')
%  xlabel('x-axis')
 
% xlabel('x-axis')
% legend('Exact soulution U ','Numerical soulution u','Initial  U0')
% title( {['\gamma =  ', num2str(gamma)]} )
% zp = BaseZoom();
% zp.plot;







% figure
%  format rat
% plot( tt,  Mass )
% xlabel('Time')
% %ylabel('Mass')
% %  legend('Mass')
% %   xticks( T0:(T-T0)/4:T )
% %  ylim([Mass(1)-1e-10  Mass(1)+1e-10])
% %  %yticks( [-1e-10,-0.5*1e-10, 0, 0.5*1e-10, 1e-10 ] )

% figure 
% plot( X , min(u,0) )
% legend('min(0,u)')
% xlabel('x')
% %ylabel('u')

% figure
% M_1 = Mass - Mass(1) ;
% plot(  tt,  M_1 );
% % ylabel('Mass-Mass_0') %已达到机器误差
% % xlabel('Time')
% % xlim([T0  T])
% % xticks( T0 : ( T- T0 )/4 :T  )
% % % ylim([-0.12* 1e-10  0.12*1e-10])
% % % yticks( [-0.12*1e-10,-0.5*0.12*1e-10, 0, 0.5*0.12*1e-10, 0.12*1e-10 ] )
% % % ylim([-1e-10  1e-10])
% % %  yticks( [-1e-10,-0.5*1e-10, 0, 0.5*1e-10, 1e-10 ] )
% % legend('Mass - Mass_0')




Emax=max(max(abs(U-u)));
E2=sqrt(sum(dx*((U-u).^2)));
format short e
% [Emax, E2] ;
fprintf('Nx = %d, Nt = %d, err_inf = %e\n',nx, nt, Emax) 
% M_err=max(max(U0))-max(max(u))  %两个结果均为正，则说明 用了限制器 结果会保界
% m_err=min(min(u))-min(min(U0))

