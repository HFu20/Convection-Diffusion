% %--------------------------------
%% 为测时间阶
function  [ m_err,M_err,u_n] = fun_FP_t(  nt,dt,T0, T,gamma,K,xmin,xmax  )  
%fun_FP_t( N_t(1), D_t(1), T0, T,gamma,K,xmin,xmax   ) ;
%  df = 10;
   dx =0.025 ;
% dt = 1e-4;
% dt = K*dx/(3*df);
nx = floor( (xmax-xmin)/dx );
T = T0 + dt*nt ;
X=( xmin:dx:xmax)';
x=X(1:end-1); 
tt = T0 : dt : T ;
Mass = zeros( length( tt ),1 ); % 初始质量也要计算
Iter = zeros( length( tt ) -1 ,1 );  % 一共nt次


%initial and Exact condition
 U0 = Ini_u( x ) ;
u_n = U0  ;
lam_n = zeros( length(u_n), 1);
kxi = 0; cx =1 ;
a = min(U0) ;
b = max(U0) ;
gprime = @(x) (-2*x +a +b) ;
 g = @(x) (b-x).*(x-a) ;

for j  = 2 : length(U0)-1
    Mass(1) = Mass(1) + U0(j)*dx ;
end

  [ matrix_Hl , matrix_Hr , matrix_Delta_x_2 , matrix_Delta_x, matrix_A12,matrix_A6] = Matrix_periodic(nx, gamma, dx, dt);
tic
for k = 1:nt

    u_nPlusHf = matrix_Hl\(matrix_Hr * u_n) + 0.5*dt*(kxi + lam_n.*gprime(u_n));  % 扩散

    [ u_nPlus1_tilde] = fun_Mstep_convec( u_nPlusHf, dt, matrix_Delta_x,matrix_A6, kxi,T0, cx,x,gamma ,K,k) ; %分成m个子步
 
    u_nPlus1 =  matrix_Hl\(matrix_Hr * u_nPlus1_tilde) + 0.5*dt*(kxi + lam_n.*gprime(u_n));  
    
         % Lagrange乘子 preserving bound and mass
%      [ kxi_plus, yinta , Iter(k) ] = fun_kxi_plus( dt, u_nPlus1,U0, kxi, dx ) ;
[ kxi_plus,count ] = lam_fun_kxi_plus( dt,  u_nPlus1,u_n,U0, kxi, dx,  lam_n  ) ;
     z = u_nPlus1 + dt*( kxi_plus - kxi - lam_n.* gprime(u_n) )  ;

    for i = 1: length( u_nPlus1 )
        if z(i) <= min(U0)
            u_nPlus1 (i) = min(U0);
        elseif z(i) >= max(U0)
            u_nPlus1 (i) = max(U0);
        else
            u_nPlus1 (i) = z(i) ;
        end
        Mass(k+1) = Mass(k+1) + u_nPlus1 (i)*dx ;
    end
   kxi = kxi_plus ;

    u_n = u_nPlus1;


end

 
% figure 
% hold on 
% plot(X,u,'o')
% plot(X,U0)
%  legend('Numeical soulution','Initial soulution')
% 
%  figure
% M_1 = Mass - Mass(1) ;
% plot( tt,  M_1 )
% ylabel('Mass- Mass_0') %已达到机器误差
% % ylim([-1e-10  1e-10])
% % yticks( [-1e-10,-0.5*1e-10,0,0.5*1e-10,1e-10 ] )
% xticks( xmin : (xmax - xmin)/4 : xmax  )
% 
%  figure
% plot( tt,  Mass )
% xlabel('Time')
% title('Mass')
%  figure
% plot( X, min(u,0) )
 % subplot(1,4,1)
% plot(x,u,'b');
% title("最终时刻的数值解")
% subplot(1,4,2)
% plot(x,U,'r');
% title("最终时刻的精确解")
% subplot(1,4,3)
% plot(x,U-u,'*-');
% subplot(1,4,4)
% % plot(x,U,'b');
% plot(x,U0,'k')
% hold on
% plot(x,u,'ro')
% legend('Exact soulution U ','Neumeical soulution u')
%     title( {['ax =   ', num2str(ax),  ' \Deltax =  ',num2str(dx),  ...
%           '  \Deltat = ', num2str(dt)] } )
% figure\
%  plot(x,u,'r-.')
% hold on
% plot(x,U0,'b')
% %legend('Initial  U0')
%    legend('数值解u','初始时刻解 U0')
%     title( {['ax =   ', num2str(ax),  ' nx =   ',num2str(nx),  ...
%           ' \Deltat =  ', num2str(dt)] } )
% 

% Emax=max(max(abs(U-u)));
% E2=sqrt(sum(dx*((U-u).^2)));
% surf(uu(1:10:end,1:25:end))
% figure
% surf(uu(1:100:end,1:25:end))
M_err=max(max(U0))-max(max(u_n)) ; %两个结果均为正，则说明 用了限制器 结果会保界
m_err=min(min(u_n))-min(min(U0));

% format short e
% [Emax, E2]
end