% clear
% %--------------------------------
% 与沈捷老师文章中给出的二阶BDF2比较
%--------------------

clear
 
%nx = 250; gamma = 0.0005;cx =1 ;
nx =600;  gamma =5e-4;  
%  nx = 1200; gamma = 2e-4;cx =1 ;

xmin = 0;
xmax = 3;
dx=((xmax  - xmin)/nx);
X=( xmin:dx:xmax)';
x=X(2:end-1); % homogeous boundary

df = 0.25;
%T0=1;T=13;
T0=1;T=4; 
% nt =80 ;
% dt =  (T-T0)/nt ;
M =1 ; dt = 1.58* dx/( 5*  df); %BDF2下 再小就不收敛了
% M =10 ; dt =   dx/( df);

% dt = 1e-5;
%    dt = dx^2;
 
nt=floor((T-T0)/dt);
tt=(T0:dt:T0+dt*nt);
T = T0+dt*nt ;

%initial and Exact condition
[U0] = Exact_exp(tt(1),  X,  gamma) ;
[U] =  Exact_exp(tt(end) , X, gamma) ;
u_minus=U0( 2 : end-1 ) ;
 


[ matrix_Hl , matrix_Hr , matrix_Delta_x_2 , matrix_Delta_x, matrix_A12,matrix_A6] = Matrix_nonp(nx-1, gamma, dx, dt) ;

%  [ matrix_Hl , matrix_Hr , matrix_Delta_x_2 , matrix_Delta_x, matrix_A1,matrix_A2] = Matrix_interi(nx-1, gamma, dx, dt) ;
Ix = eye( nx, nx ) ;

a = min(U0) ;
b = max(U0) ;
lam_n =zeros(length(u_minus),1) ;
gprime = @(x) (-2*x +a +b) ;
 g = @(x) (b-x).*(x-a) ;
lam = 0;

kxi = 0;
Mass_u = zeros( length( tt ),1 );
for j  = 1 : length(U0) 
    Mass_u(1) = Mass_u(1) + U0(j)*dx ;
end
%BDF1
T_lag = tic ;
k = 1 ;% t = t0+dt
Left_1 = matrix_A12*matrix_A6 - dt*gamma*matrix_A6*matrix_Delta_x_2 ; 
u_0 =  Left_1\( matrix_A12*matrix_A6*u_minus - dt*matrix_A12*matrix_Delta_x*fun_f(u_minus)  ) ;
%Mass_u(k+1) =sum(u_0 )*dx ;

%BDF2
Left_2 = matrix_A12*matrix_A6* 3/(2 )  - gamma*dt*matrix_A6*matrix_Delta_x_2 ;
for k = 2:nt

u_plus =  Left_2\( 2*matrix_A12*matrix_A6* u_0 - 1/2*matrix_A12*matrix_A6* u_minus ...
         - dt* matrix_A12*matrix_Delta_x*fun_f(2*u_0-u_minus) ); %  ...
   %     + dt*matrix_A12*matrix_A6* lam*gprime(u_0) ) ;  %     %    
 %    +  dt*matrix_A12*matrix_A6*(kxi + lam_n.*gprime(u_0) ) ); 
% yinta = - dt*2/3*lam*gprime(u_0) ;
% z = yinta + u_plus ;
% if z<= a
%     u = a ;
%      lam = 3*( a - (u_plus + yinta) )/( 2*dt*gprime(a) ) ;
% elseif z>=b
%     u = b ;
%      lam = 3*( b - (u_plus + yinta) )/( 2*dt*gprime(b) ) ;
% else 
%     u = z ;
%     lam = 0;
% end
% u_plus = u ;

%          % Lagrange乘子 preserving bound and mass
% [ kxi_plus(k),count(k) ] = lam_fun_kxi_plus( dt,  u_plus ,u_0,U0, kxi, dx,  lam_n  ) ;
% yinta =  1/1*dt*(  kxi_plus(k) - kxi -lam_n.*g(u_0) ) ;    
%      z = u_plus  + yinta  ;
% kxi = kxi_plus(k) ;
%     for i = 1: length( u_plus  )
%         if z(i) <= min(U0)
%             u_plus  (i) = min(U0);
%               lam_n(i) = 1*(a - (u_plus (i) + yinta(i)))/(dt*gprime(a))  ;
%         elseif z(i) >= max(U0)
%             u_plus  (i) = max(U0);
%              lam_n(i) =1* (b - (u_plus (i) + yinta(i)))/(dt*gprime(b))  ;
%         else
%             u_plus  (i) = z(i) ;
%                lam_n(i) = 0 ;
%         end
%     end
% Mass_u(k+1) =sum(u_plus )*dx ;

    u_minus = u_0;
    u_0 = u_plus ;

end
     for i = 1: length( u_plus ) 
        if u_plus(i) <= min(U0)
            u_plus (i) = min(U0);  
         elseif u_plus (i) >= max(U0)
            u_plus (i) = max(U0);
        end  
     end

u = U ;
u( 2:end-1 ) = u_plus;
 toc(T_lag) ;
 min(u)
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
plot( X ,u,'b-d')
legend('Exact soulution u ','Numerical soulution U')
 xlabel('x-axis')

figure
plot( X ,U-u,'b--','LineWidth',1.3);
ylim([ -6e-3 6e-3 ] ) 

%  figure
%  Mass_err_u = Mass_u -Mass_u(1) ;
% plot(tt,Mass_err_u)
%  xlabel('Time')
% Emax=max(max(abs(U-u)));
% E2=sqrt(sum(dx*((U-u).^2)));
% surf(uu(1:10:end,1:25:end))
% figure
% surf(uu(1:100:end,1:25:end))
% M_err=max(max(U0))-max(max(u_minus))  %两个结果均为正，则说明 用了限制器 结果会保界
% m_err=min(min(u_minus))-min(min(U0))

% format short e
%  [Emax, E2]