clear

a=1.5;c=0.5;%精确解所需参数（用于边界）
xspan = [-1, 3]; tspan = [0, 1];
T =  tspan(2)-tspan(1);
% gamma=1/64; % 对应S = 0
gamma = 1e-3 ; % 对应S非0
df=2;
%==================================
%Conv + Diffu
%  dx =gamma/(2*df);miu=1/3;dt=gamma/(3*df); nx= floor(( xspan(2)  - xspan(1))/dx);
%--------------------
%Only convection
 nx=1000; dx=(( xspan(2)  -  xspan(1))/nx);  dt = dx/( 3*df);
%原来   nx = 200; nt = 100;
%========================================
nt =floor((tspan(2) - tspan(1)) / dt);
%     dx = abs(xspan(2) - xspan(1)) / (nx + 1); % x方向的步长
%     dt = abs(tspan(2) - tspan(1)) / nt; % t方向的步长
% 计算网格节点
x_all = xspan(1):dx: xspan(1)+nx*dx ;  % x方向所有节点
t_all = tspan(1):dt: tspan(1)+nt*dt ; % t方向所有节点
x_in = x_all(2:end -1)'; % x方向内点
t_in = t_all(1:end -1); % t方向内点
%     [x, t] = meshgrid(x_in, t_in); % 网格节点
b=gamma;
% u(x,t)=a-c*tanh(c/(2*b)*(x-a*t));
% initial and Exact sol.
U0 = Exact_sol( x_all' , tspan(1) ,gamma,a,c ) ;
U   = Exact_sol( x_all' , tspan(2),gamma,a,c ) ;
u_n=U0( 2: end-1 );

[matrix_LHS_step1, matrix_RHS_step1, matrix_A6, matrix_Delta_x] = Matrix_nonp (nx-1 ,dt, dx, gamma); % generate nx+1 * nx+1 matrix

G_n=zeros(nx-1,1);
G_nplus1=zeros(nx-1,1);
% G_deltax=zeros(nx+2,1);
% G_deltax(1)=1;  % From delta_x
% G_deltax(end)=-1;
g1 = zeros(nx-1,1);
g2 =   zeros(nx-1,1);
for k = 1:nt
    %boundary
    G_n(1) =  Exact_sol( x_all(1), t_all(k), gamma, a, c ) ; % 第k层的u0;
    G_n(end)= Exact_sol( x_all(end), t_all(k), gamma, a, c ) ; % 第k层的uend

    G_nplus1(1) =  Exact_sol( x_all(1), t_all(k+1), gamma, a, c ) ;%第k+1层的u0;
    G_nplus1(end)= Exact_sol( x_all(end), t_all(k+1), gamma, a, c ); % 第k+1层的uend;

    u_nPlusHf = matrix_LHS_step1\( (matrix_RHS_step1 * u_n)+...
        (matrix_RHS_step1(2,1)- matrix_LHS_step1(2,1))*G_n ) ;
    %这是由于三对角+对称 matrix_RHS_step1(2,1) = matrix_RHS_step1(1,2),LHS_step1也是

     g1(1) = - matrix_A6(2,1)*G_nplus1(1) + matrix_A6(2,1)*G_n(1)- dt* matrix_Delta_x(2,1)*fun_f( G_n(1) ) + dt*matrix_A6(2,1)*fun_Source_wan(x_all(1), t_all( k )) ;
    g1(end) = -matrix_A6(1,2)*G_nplus1(end) + matrix_A6(2,1)*G_n(end)- dt* matrix_Delta_x(1,2)*fun_f( G_n(end) ) + dt*matrix_A6(1,2)*fun_Source_wan(x_all(end), t_all( k )) ;
    u_nPlusHfStar = matrix_A6\( matrix_A6*u_nPlusHf   - dt * matrix_Delta_x  * fun_f( u_nPlusHf  ) + dt*matrix_A6*fun_Source_wan( x_in,t_all(k)  ) +  g1  );

    g2(1) = -matrix_A6(2,1)*G_nplus1(1) + matrix_A6(2,1)*G_n(1) - dt/2* matrix_Delta_x(2,1)* ( fun_f( G_n(1) ) + fun_f( G_nplus1(1) ) )...
        + dt/2*matrix_A6(2,1)*( fun_Source_wan(x_all(1), t_all( k ))+ fun_Source_wan(x_all(1), t_all( k +1))  );
    g2(end) = -matrix_A6(1,2)*G_nplus1(end) + matrix_A6(1,2)*G_n(end)- dt/2* matrix_Delta_x(1,2) * ( fun_f( G_n(end) ) + fun_f( G_nplus1(end) ))...
        + dt/2*matrix_A6(1,2)* ( fun_Source_wan(x_all(end), t_all( k )) + fun_Source_wan(x_all(end), t_all( k+1 )) ) ;
    u_nPlus1_tilde = matrix_A6\(  matrix_A6*u_nPlusHf  - 0.5 *dt * matrix_Delta_x * ( fun_f(u_nPlusHfStar)+ fun_f(u_nPlusHf) ) ...
        + dt/2*matrix_A6*( fun_Source_wan( x_in,t_all(k) )  + fun_Source_wan( x_in,t_all(k+1) ) )+ g2 );

    u_nPlus1 = matrix_LHS_step1\( (matrix_RHS_step1 * u_nPlus1_tilde)+...
        matrix_RHS_step1(2,1)*G_n- matrix_LHS_step1(2,1)*G_nplus1);

    u_n = u_nPlus1;

end
     for i = 1: length( u_nPlus1 ) 
        if u_nPlus1(i) <= min(U0)
            u_nPlus1 (i) = min(U0);  
         elseif u_nPlus1 (i) >= max(U0)
            u_nPlus1 (i) = max(U0);
        end  
     end
      u_n = u_nPlus1; 
u = U ;
u(2:end-1) = u_n;



% figure
% plot(x_all,U,'b');
% hold on
% plot(x_all,u,'*')
% plot(x_all,U0)
% title( {['\gamma =  ', num2str(gamma)]} )
% legend('Exact solution','Numerical solution','Initial condition')
 
%
figure
% ax1 = axes('Position', [0.1 0.1 0.8 0.8]);
%  axes(ax1) %大窗口绘图
plot(x_all,U,'b-');
hold on 
plot(x_all,u,'*')
%plot(x_all,U0)
xlabel('x-axis')
title( '\gamma = 1e-3'  )
legend('Exact solution','Numerical solution')
%   magnify
% ylim( [0.8 2.2])
%  
% ax2 = axes('Position',[0.5 0.34 0.3 0.3]);
% axes(ax2) %小窗口绘图
%   a = 610 : 618; % nx = 1000
% %a = 621: 628; %上拐角
% % a = 641: 660;  %下拐角
% %a = 730 : 775; % nx =1200
%  x1 = x_all(a) ;
% plot(x1,U(a),'b-');
% hold on 
% plot(x1,u(a),'*')
% u1 = min( u(a) ) ;u2 = max( u(a) ) ;
% xlim( [x1(1) x1(end)])
% ylim( [u1, u2])

Emax=max(max(abs(U-u)));
E2=sqrt(sum(dx*((U-u).^2)));

fprintf('lim:nx =%d,gamma=%f,nt=%f\n',nx,gamma,nt)
M_err = max(U0)-max(u) %两个结果均为正，则说明 用限制器以后 结果保界
m_err = min(u)-min(U0)
fprintf('[Emax, E2] =[ %e, %e ]\n',[Emax, E2])

