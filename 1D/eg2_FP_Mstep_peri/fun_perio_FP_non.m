% %--------------------------------
% 只给初值 U0= sin(2*pi*x)+1;
% 针对1D非线性burgers方程
%------------------
function [u_n ,Mass_1, tt, max_u, min_u,count, Energy] = fun_perio_FP_non(  x ,  nx, dx, nt,T ,T0, gamma )

dt = ( ( T - T0 )/nt )  ; 
%initial  condition
 U0 = Ini_u( x ) ;
u_n = U0;

[ matrix_Hl , matrix_Hr , matrix_Delta_x_2 , matrix_Delta_x, matrix_A1,matrix_A2] = Matrix_periodic(nx, gamma, dx, dt) ;% generate nx+1 * nx+1 matrix

kxi = 0;
T = T0 + nt*dt ;
tt = T0 : dt : T ;
count = 0 ;

Mass_u = zeros( length( tt ),1 ); % 初始质量也要计算
Count = zeros( length( tt ) -1 ,1 );  % 一共nt次
Energy = zeros( length( tt ),1 );

%Energy(1) =   fun_Energy( x, U0, dx  ) ;
Mass_u(1) = sum(U0)*dx ;

M = 1; cx =1;
kxi_plus = zeros(nt,1) ;
lam_n = 0 ;
a= min(U0) ;
b = min(U0) ;
gprime = @(x) (-2*x +a +b) ;
 g = @(x) (b-x).*(x-a) ;

tic
for k = 1:nt

    u_nPlusHf = matrix_Hl\(matrix_Hr * u_n) + 0.5*dt*(kxi + lam_n.*gprime(u_n));   % 扩散

   [ u_nPlus1_tilde] = fun_Mstep_convec( u_nPlusHf, dt, matrix_Delta_x,matrix_A2, kxi,T0, cx,x,gamma ,M,k) ; %分成m个子步

    u_nPlus1 = matrix_Hl\(matrix_Hr * u_nPlus1_tilde) + 0.5*dt*(kxi + lam_n.*gprime(u_n));   % 扩散
    
%          % Lagrange乘子 preserving bound and mass
% [ kxi_plus(k),count ] = lam_fun_kxi_plus( dt,  u_nPlus1,u_n,U0, kxi, dx,  lam_n  ) ;
% yinta =  1/2*dt*(  kxi_plus(k) - kxi -lam_n.*g(u_n) ) ;    
%      z = u_nPlus1 + yinta  ;
% kxi = kxi_plus(k) ;
%     for i = 1: length( u_nPlus1 )
%         if z(i) <= min(U0)
%             u_nPlus1 (i) = min(U0);
%         elseif z(i) >= max(U0)
%             u_nPlus1 (i) = max(U0);
%         else
%             u_nPlus1 (i) = z(i) ;
%         end
%     end

 Mass_u(k+1) =sum(u_nPlus1)*dx ;
  Energy(k+1) = fun_Energy( x, u_nPlus1,dx  ) ;
 u_n = u_nPlus1;


end
toc

% M_err = max(max(U0))-max(max(u_n))  %两个结果均为正，则说明 用了限制器 结果会保界
% m_err = min(min(u_n))-min(min(U0)) 
Mass_1 = Mass_u - Mass_u(1) ;

max_u = max(u_n);
min_u = min( u_n ) ;
end

% format short e
% [Emax, E2]