% %--------------------------------
% 只给初值 U0= sin(2*pi*x)+1;
% 针对1D非线性burgers方程
%------------------
function [u ,U0, Mass_1, tt, max_u, min_u] = fun_perio_U0_sin1( x ,  nx_1, dx, dt,T ,T0, gamma )

nt =floor( ( T - T0 )/dt )  ; 
%initial  condition
 U0 = Sin_sol( x' );
u_n=U0;

[matrix_Hl,matrix_Hr, matrix_A2, matrix_Delta_x] =  Matrix(nx_1 ,dt, dx, gamma); % generate nx+1 * nx+1 matrix

kxi = 0;
tt = T0 : dt : T ;
Mass = zeros( length( tt ),1 ); % 初始质量也要计算
Count = zeros( length( tt ) -1 ,1 );  % 一共nt次

for j  = 1 : length(U0)
    Mass(1) = Mass(1) + U0(j)*dx ;
end

M = 10; cx =1;
tic
for k = 1:nt

    u_nPlusHf = matrix_Hl\(matrix_Hr * u_n);   % 扩散

         [ u_nPlus1_tilde] = fun_Mstep_convec( u_nPlusHf, dt, matrix_Delta_x,matrix_A2, kxi,T0, cx,x,gamma ,M,k) ; %分成m个子步
%% 未加
%   u_nPlusHfStar =u_nPlusHf - matrix_A2\( dt * matrix_Delta_x  * fun_f(u_nPlusHf)  ) + ...
%         dt*fun_Source(x, T0+(k-1)*dt, gamma,cx )+dt*kxi ;  %对流
% 
%     u_nPlus1_tilde = u_nPlusHf - matrix_A2\(  0.5 *dt * matrix_Delta_x  ...
%         * (fun_f(u_nPlusHf) + fun_f(u_nPlusHfStar))  )   +  ...
%         dt/2* ( fun_Source(x, T0+ (k-1)*dt, gamma,cx ) + fun_Source( x, T0+(k)*dt, gamma,cx ) )+ dt*kxi ;   %对流

   %-------------------------limter--------------------------------
%     u_nPlus1_tilde  =limiter_2sides_final(nx,u_nPlus1_tilde,min(U0),max(U0));
%-----------------------------------------------------------------------------

    u_nPlus1 = matrix_Hl\(matrix_Hr * u_nPlus1_tilde);
    
             % Lagrange乘子 preserving bound and mass
     [ kxi_plus, yinta , Count(k) ] = fun_kxi_plus( dt, u_nPlus1,u_n, U0, kxi, dx ) ;
     z = u_nPlus1 + yinta  ;

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

u = u_n;

M_err = max(max(U0))-max(max(u))  %两个结果均为正，则说明 用了限制器 结果会保界
m_err = min(min(u))-min(min(U0)) 
% Mass_1 = Mass - Mass(1) ;
 Mass_1 = Mass ;
max_u = max(u_n);
min_u = min( u_n ) ;
end

% format short e
% [Emax, E2]