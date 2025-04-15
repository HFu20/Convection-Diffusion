%% 23. 03.24 
%求空间阶数 但精度不理想
% 取 dt = miu* dx^2 同时验证了时间精度

%%================================================
% Eq: u_t+ cx*u_x=gamma *u_{xx}
%+++++++++++++++++++++++++++++++++++++++++++
%% Exact solution:  u(x,y,t)=1/sqrt(1+t)*exp(-(x-a*(1+t)).^2/(4*b*(1+t));
clear
% %------------参数----------------
%%%%  
  gamma=0.005;  

%     N = [  512 1024 2056 ];
N = [ 180  360 540  720 ];
 xmin =0;xmax = 3; %gauss
% xmin =0;xmax = 2;
 T0=1;T_old =4;
 
%===================================

nstep_x = length(N);

Error_Table = zeros( nstep_x ,5 );
Bound_Table =  zeros( nstep_x ,3 );
    

for i = 1:nstep_x
    
  nx = N(i);
 M = 10    ;
  dx= (xmax  - xmin)/nx; 
    dt = 1e-4; 
%  dt = dx^2 ;
  nt=floor((T_old-T0)/dt);
% nt = 1000;
% % dt =  (T-T0)/nt;
%  dt =dx/(1*df);
%   nt=floor((T-T0)/dt);
X=( xmin:dx:xmax)';
x=X(2:end-1); % homogeous boundary
 

t=(T0:dt:T0+dt*nt);
T = T0+dt*nt ;
if T ~= T_old
    fprintf(1,'\n [T, T_old] = [ %f, %f]\n',T , T0  );
end
 %  [ E2,Emax,m_err,M_err,u_n] = fun_eg3_exp_diff_homo( dx,dt, nx, nt, T0, T, X,x,gamma ) ;
[E2,Emax,m_err,M_err,u_n,time(i)] = lam_fun_exp_diff_homo( dx,dt, nx, nt, T0, T, X,x,gamma , M) ; %额外使用了lamda的信息

    % put n in the 1st column of table
    Error_Table(i,  1) = nx ;
    % compute the Lmax norm of error and put it in the 2nd column
    Error_Table(i, 2)  =  Emax ;
    % compute the L2 norm of error and put it in the 2nd column
    Error_Table(i, 4)  =  E2 ;

    Bound_Table(i,  1) = nx ;
    Bound_Table(i, 2) =  M_err;
    Bound_Table( i, 3 ) =  m_err;

    if i ~=1
        Error_Table(i, 3) = log(Error_Table(i, 2)/Error_Table(i-1, 2))...
            /-log(Error_Table(i, 1)/Error_Table(i-1, 1));

        % compute the ratio for E2, and put it in the 5th column of table
        Error_Table(i, 5) = log(Error_Table(i, 4)/Error_Table(i-1, 4))...
            /-log(Error_Table(i,1)/Error_Table(i-1, 1));

    end


end

format short e
% convert array to table for visualization
TableResults_Err = array2table(Error_Table, 'VariableNames', ...
    {'Nx', 'E_max', 'ratio E_max', 'E_2', 'ratio E_2'})
%--------------------------------------------------
TableResults_Bound = array2table(Bound_Table, 'VariableNames', {'Nx', 'M_err', 'm_err'})

writetable(Error_Table, 'TableResults_Err.xls')
writetable(Bound_Table, 'TableResults_Bound.xls')


%  ax1 = axes('Position', [0.1 0.1 0.8 0.8]);
%   axes(ax1) %大窗口绘图
[U] =  Exact_exp(t(end) , X, gamma) ;
% figure
% plot( X ,u_n,'ro')
% hold on 
%  plot( X  ,U,'k')
% xlabel('x-axis')
% legend( 'Numerical soulution u')
% fprintf(1,'\n time = [ %f, %f,%f,%f,%f,%f]\n',time  );
 
