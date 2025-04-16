% %--------------------------------
% 只给初值 U0= sin(2*pi*x)+1;
% 针对1D非线性burgers方程
%------------------
function [u_n ,U0, Mass_1, tt, max_u, min_u] = fun_perio_FP_U0_sin1(  x ,  nx, dx, nt,T ,T0, gamma )

dt = ( ( T - T0 )/nt )  ; 
%initial  condition
 U0 = Ini_u( x ) ;
u_n = U0;

[ matrix_Hl , matrix_Hr , matrix_Delta_x_2 , matrix_Delta_x, matrix_A1,matrix_A2] = Matrix_periodic(nx, gamma, dx, dt) ;% generate nx+1 * nx+1 matrix

kxi = 0;
T = T0 + nt*dt ;
tt = T0 : dt : T ;


Mass = zeros( length( tt ),1 ); % 初始质量也要计算
Count = zeros( length( tt ) -1 ,1 );  % 一共nt次

for j  = 1 : length(U0)
    Mass(1) = Mass(1) + U0(j)*dx ;
end

M = 5; cx =1;
tic
for k = 1:nt

    u_nPlusHf = matrix_Hl\(matrix_Hr * u_n);   % 扩散

   [ u_nPlus1_tilde] = fun_Mstep_convec( u_nPlusHf, dt, matrix_Delta_x,matrix_A2, kxi,T0, cx,x,gamma ,M,k) ; %分成m个子步

    u_nPlus1 = matrix_Hl\(matrix_Hr * u_nPlus1_tilde);
    
             % Lagrange乘子 preserving bound and mass
%      [ kxi, yinta , Count(k) ] = fun_kxi_plus( dt, u_nPlus1,U0, kxi, dx ) ;
%      z = u_nPlus1 + yinta  ;

     for i = 1: length( u_nPlus1 )
%         if z(i) <= min(U0)
%             u_nPlus1 (i) = min(U0);
%         elseif z(i) >= max(U0)
%             u_nPlus1 (i) = max(U0);
%         else
%             u_nPlus1 (i) = z(i) ;
%         end
        Mass(k+1) = Mass(k+1) + u_nPlus1 (i)*dx ;
    end
 

    u_n = u_nPlus1;

end

% M_err = max(max(U0))-max(max(u_n))  %两个结果均为正，则说明 用了限制器 结果会保界
% m_err = min(min(u_n))-min(min(U0)) 
Mass_1 = Mass - Mass(1) ;
%  Mass_1 = Mass ;
max_u = max(u_n);
min_u = min( u_n ) ;
end

% format short e
% [Emax, E2]