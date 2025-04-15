% %--------------------------------
%exact solution: 
% (x/t).*exp(-(x.*x-0.25*t)/(4*ax*t))./(sqrt(t)+exp(-(x.*x-0.25*t)/(4*ax*t)));
%--------------------------
%% 2023.03.20 
%不能用周期 
% 额外用了lamda
%===================================
function  [E2,Emax,m_err,M_err,u,t] = lam_fun_exp_diff_homo( dx,dt, nx, nt, T0, T, X,x,gamma, M )
cx =1 ;
%initial and Exact condition
[U0] = Exact_exp(T0,  X,  gamma) ;
[U] =  Exact_exp(T , X, gamma) ;
u_n=U0( 2 : end-1 ) ;
kxi = 0;
lam_n = zeros( length(u_n), 1);
tt = T0 : dt : T ;
Mass = zeros( length( tt ),1 ); % 初始质量也要计算
Iter = zeros( length( tt ) -1 ,1 );  % 一共nt次

a = min(U0) ;
b = max(U0) ;
gprime = @(x) (-2*x +a +b) ;
for j  = 2 : length(U0)-1
    Mass(1) = Mass(1) + U0(j)*dx ;
end

[ matrix_Hl , matrix_Hr , matrix_Delta_x_2 , matrix_Delta_x, matrix_A1,matrix_A2] = Matrix_nonp(nx-1, gamma, dx, dt) ;

% M = 1 ;
Time = tic ;
for k = 1:nt

    u_nPlusHf = matrix_Hl\(matrix_Hr * u_n)  + dt/2*(kxi + lam_n.*gprime(u_n)) ;   % 扩散    

%     [ u_nPlus1_tilde] =...
%         fun_Mstep_convec( u_nPlusHf, dt, matrix_Delta_x,matrix_A2, kxi,T0, cx,x,gamma ,M,k,lam_n,gprime,u_n) ; %分成m个子步
  [ u_nPlus1_tilde] =...
         fun_Mstep_convec( u_nPlusHf, dt, matrix_Delta_x,matrix_A2, kxi,T0, cx,x,gamma ,M,k) ; %分成m个子步

    u_nPlus1 =matrix_Hl\(matrix_Hr * u_nPlus1_tilde)  + dt/2*(kxi + lam_n.*gprime(u_n)) ;   % 扩散

      % %--------------------------------------------
      %  Lagrange乘子 preserving bound and mass
%       [ kxi, yinta , Iter(k) ] = fun_kxi_plus( dt, u_nPlus1, u_n,U0, kxi, dx ) ;
%      z = u_nPlus1 + yinta  ;

% [ kxi, yinta,count ] = lam_fun_kxi_plus( dt,  u_nPlus1,u_n,U0, kxi, dx, lam_n ) ;
% z= u_nPlus1 + yinta  ;
%            for i = 1: length( u_nPlus1 )
% % % % % 
% % % % %         if z(i) <= min(U0)
% % % % %             u_nPlus1 (i) = min(U0);  
% % % % %             lam_n(i) =2* (a-z(i) )/(dt*gprime(a)) ;
% % % % %         elseif z(i) >= max(U0)
% % % % %             u_nPlus1 (i) = max(U0);
% % % % %              lam_n(i) = 2*(b-z(i) )/(dt*gprime(b)) ;
% % % % %         else
% % % % %             u_nPlus1 (i) = z(i) ;
% % % % %              lam_n(i) =0 ;
% % % % %         end     
% % % % 
%         if u_nPlus1(i) <= min(U0)
%             u_nPlus1 (i) = min(U0);  
%          elseif u_nPlus1 (i) >= max(U0)
%             u_nPlus1 (i) = max(U0);
%         end
% % % %       Mass(k+1) = Mass(k+1) + u_nPlus1 (i)*dx ;

%          end

    %--------------------------------------------

    u_n = u_nPlus1;

end
%      for i = 1: length( u_nPlus1 ) 
%         if u_nPlus1(i) <= min(U0)
%             u_nPlus1 (i) = min(U0);  
%          elseif u_nPlus1 (i) >= max(U0)
%             u_nPlus1 (i) = max(U0);
%         end  
%      end
%       u_n = u_nPlus1;

u = U ;
u( 2:end-1 ) = u_n;
t = toc(Time)  ;


% E1= sum(abs(U-u)*dx) ;
Emax=max(max(abs(U-u)));
E2=sqrt(sum(dx*((U-u).^2)));
format short e
% [E1 Emax, E2];

M_err=max(max(U0))-max(max(u));  %两个结果均为正，则说明 用了限制器 结果会保界
m_err=min(min(u))-min(min(U0));

end
