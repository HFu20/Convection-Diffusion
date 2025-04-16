% %--------------------------------
%% 23.03.24 保界保质量
clear
%  nx=500;
% nt=1000;
% % N = [40 50 60 70 80];
 gamma = 1  ;
df = 1 ;

xmin = - 2*pi;
xmax = 2*pi;
T0=0;
T=2;
 
%  dx = 0.02 ;
% nx=floor ((xmax  - xmin)/dx);
nt =40 ;  
dt=((T-T0)/nt);

nx =81;
dx= ((xmax  - xmin)/nx);
% nt =100 ; %为了能取到最终时刻
% dt=((T-T0)/nt);
%  dt=dx/(6*df);
% nt = floor( (T-T0)/dt );
T =T0+nt*dt;
t=(T0:dt:T);

X=( xmin:dx:xmax)';
x=X(1:end-1);


 U0 = Ini_u( x ) ;
u_n = U0  ;


[ matrix_Hl , matrix_Hr , matrix_Delta_x_2 , matrix_Delta_x, matrix_A12,matrix_A6] = Matrix_periodic(nx, gamma, dx, dt) ;
%  [ matrix_Hl , matrix_Hr , matrix_Delta_x_2 , matrix_Delta_x, matrix_A1,matrix_A2] = Matrix_interi(nx-1, gamma, dx, dt) ;


kxi = 0;
tt = T0 : dt : T ;
Mass_u = zeros( length( tt ),1 ); % 初始质量也要计算
Mass_U = zeros( length( tt ),1 ); % 初始质量也要计算
 Energy = zeros( length( tt ),1 );
count = zeros( nt,1 );


Mass_u(1)  = sum(U0)*dx ;
 Energy(1) =   fun_Energy( x, U0, dx  ) ;

M =1; cx =1;
kxi_plus = zeros(nt,1) ;
lam_n =zeros(length(u_n),1) ;
a= min(U0) ;
b = max(U0) ;
gprime = @(x) (-2*x +a +b) ;
 g = @(x) (b-x).*(x-a) ;
tic
for k = 1:nt

    u_nPlusHf = matrix_Hl\(matrix_Hr * u_n + 0.5*matrix_A12*matrix_A6* dt*(kxi + lam_n.*gprime(u_n)) );   % 扩散

    [ u_nPlus1_tilde] = fun_Mstep_convec( u_nPlusHf, dt, matrix_Delta_x,matrix_A6, kxi,T0, cx,x,gamma ,M,k) ; %分成m个子步

    u_nPlus1 =  matrix_Hl\(matrix_Hr * u_nPlus1_tilde+ 0.5*matrix_A12*matrix_A6* dt*(kxi + lam_n.*gprime(u_n)) ); 

         % Lagrange乘子 preserving bound and mass
[ kxi_plus(k),count(k) ] = lam_fun_kxi_plus( dt,  u_nPlus1,u_n,U0, kxi, dx,  lam_n  ) ;
yinta =  1/1*dt*(  kxi_plus(k) - kxi -lam_n.*g(u_n) ) ;    
     z = u_nPlus1 + yinta  ;
kxi = kxi_plus(k) ;
    for i = 1: length( u_nPlus1 )
        if z(i) <= min(U0)
            u_nPlus1 (i) = min(U0);
              lam_n(i) = 1*(a - (u_nPlus1(i) + yinta(i)))/(dt*gprime(a))  ;
        elseif z(i) >= max(U0)
            u_nPlus1 (i) = max(U0);
             lam_n(i) =1* (b - (u_nPlus1(i) + yinta(i)))/(dt*gprime(b))  ;
        else
            u_nPlus1 (i) = z(i) ;
               lam_n(i) = 0 ;
        end
    end

 Mass_u(k+1) =sum(u_nPlus1)*dx ;
 Energy(k+1) = fun_Energy( x, u_nPlus1,dx  ) ;
 u_n = u_nPlus1;


end
toc
 
% writematrix(u_n,'u_n.xls')
%  readmatrix('u_n.xls')
% 
% M = readmatrix('u_n.xls');

figure 
hold on 
 plot(x,u_n,'-s')
plot(x,U0,'-o')
 legend('Numeical soulution','Initial solution' )
xlim([xmin xmax])


figure
M_1 = Mass_u - Mass_u(1) ;
plot( tt,  M_1,'r-d' )
ylabel('Mass_{err}') %已达到机器误差
 ylim([-1e-13   1e-13])
% yticks( [-1e-10,-0.5*1e-10,0,0.5*1e-10,1e-10 ] )
xlim([ 0,  T])
xlabel('Time')

figure
% E_err =  Energy   - Energy(1) ;
plot( tt,    (Energy),'-.','LineWidth',1.5 )
ylabel('Entropy') %已达到机器误差
xlabel('Time')
%  ylim([-1e-14  1e-14])
% yticks( [-1e-10,-0.5*1e-10,0,0.5*1e-10,1e-10 ] )
xlim([ 0,  T])
 

% Emax=max(max(abs(U-u)));
% E2=sqrt(sum(dx*((U-u).^2)));
% surf(uu(1:10:end,1:25:end))
% figure
% surf(uu(1:100:end,1:25:end))
% M_err=max(max(U0))-max(max(u_n))  %两个结果均为正，则说明 用了限制器 结果会保界
% m_err=min(min(u_n))-min(min(U0))

% format short e
% [Emax, E2]