function   [Mass_u,u_n, kxi_plus,yinta, Gplus, Mass_err_u] = mass_Solve_u_inter_DYak_S( miu_x, miu_y,  lam_x, lam_y, nx , ny,dx, dy,  gamma, nt, dt, X, Y,x,y, U0, T0,T,M)

%----matrix -------------
[ matrix_Delta_x_2,matrix_Delta_x ,matrix_A1x ,matrix_A2x ] = Matrix_periodic(nx) ;
[ matrix_Delta_y_2,matrix_Delta_y ,matrix_A1y ,matrix_A2y ] = Matrix_periodic(ny);

H_1x = matrix_A1x - gamma*miu_x/4*matrix_Delta_x_2;
H_1y = matrix_A1y - gamma*miu_y/4*matrix_Delta_y_2;
H_2x = matrix_A1x + 1/4*gamma*miu_x*(matrix_Delta_x_2);
H_2y = matrix_A1y + 1/4*gamma*miu_y*(matrix_Delta_y_2);

u_n = U0 ;
a = min( min(U0) );
b = max(max(U0)) ;
lam_n = zeros(length(u_n) ) ;
yinta = zeros(length(u_n) ) ;
gprime = @(x) (-2*x +a +b) ;
 g = @(x) (b-x).*(x-a) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
% kxi_plus = -1e-12 ;
kxi = 0 ;
tt = T0 : dt : T ;
Mass_u = zeros( length( tt ),1 ); % 初始质量也要计算
 Mass_S = zeros( length( tt ),1 ); % 初始质量也要计算
% Count = zeros( nt ,1);
kxi_plus = zeros( nt ,1 );
Gplus = zeros( nt ,1 );
% Gminus = zeros( nt ,1 );

%初始质量
 Mass_u(1) =  sum( sum(u_n ) )  *dx*dy ; % 只考虑内点质量（边界为0）

 for k = 0 : nt
 Mass_S(1:end) =    sum(sum( fun_Source( X, Y, gamma , T0+k*dt )  ) )*dx*dy ;
 end
% Mass_S(1) = 0;
%  yinta只跟t有关，故第k层时，只是一个数
Mass_S_half = 0.5*dt*( Mass_S(2:end) + Mass_S(1:end-1) ) ;
%%%%%%%%%%%%%%%%%%%

 
tic
for k = 1:nt

    %保留A_1y*unPlusHf1,第2步里直接运用
    %     u_bar_nPlusHf1 =  H_1x\( H_2x*u_n*matrix_A1y'+ miu_y *gamma/ 2*matrix_A1x*u_n*matrix_Delta_y_2'  );
    %
    %     u_nPlusHf2       = (u_bar_nPlusHf1-miu_y*gamma/4* u_n*matrix_Delta_y_2')/ H_1y' ;

    u_bar_nPlusHf1 =  H_1x\( H_2x*u_n*H_2y'  ) + 0.5*dt*(kxi + lam_n.*gprime(u_n)); 
    u_nPlusHf2       = u_bar_nPlusHf1/ H_1y' ;

%      [  u_nPlusHf4 ] = fun_Mstep_convec( u_nPlusHf2, matrix_A2x, matrix_Delta_x, lam_x, lam_y, matrix_Delta_y, matrix_A2y, dt,X, Y, gamma, T0,  k, kxi ,M ) ;
 [  u_nPlusHf4 ] =   lam_fun_Mstep_convec( u_nPlusHf2, matrix_A2x, matrix_Delta_x, lam_x, ...
                                  lam_y, matrix_Delta_y, matrix_A2y, dt,X, Y, gamma, T0,  k,  M);
% lam_fun_Mstep_convec( u_nPlusHf2, matrix_A2x, matrix_Delta_x, lam_x, ...
%                                  lam_y, matrix_Delta_y, matrix_A2y, dt,X, Y, gamma, T0,  k,  M ,kxi,lam_n, u_n,gprime );
    u_bar_nPlusHf5 =  H_1x\( H_2x*u_nPlusHf4*H_2y'  )  + 0.5*dt*(kxi + lam_n.*gprime(u_n)) ;
    u_nPlus1      = u_bar_nPlusHf5/ H_1y'  ; 
    %%--------------------------------------------

[ kxi_plus(k) ] = new_fun_kxi_plus( dt, u_nPlus1, u_n, a, b , kxi, dx,dy,Mass_S_half(k),lam_n ) ;
yinta =  1/1*dt*(  kxi_plus(k) - kxi - lam_n.*g(u_n) ) ;    
kxi = kxi_plus(k) ;
z = u_nPlus1+ yinta  ;
%  
    for i = 1: nx
        for j = 1: ny
            if z( i,j ) <= a
                u_nPlus1 ( i,j ) = a;
                 lam_n(i,j) =1* (a-z(i,j) )/(dt*gprime(a)) ;
            elseif z( i,j ) >= b
                u_nPlus1 (i,j) = b;
                 lam_n(i,j) = 1*(b-z(i,j) )/(dt*gprime(b)) ;
            else
                u_nPlus1 (i,j) = z(i, j ) ;
                 lam_n(i,j) =0 ;
            end
%             if u_nPlus1( i,j ) <= min( min(U0) )
%                 u_nPlus1 ( i,j ) = min(min(U0));
%             elseif u_nPlus1( i,j ) >= max(max(U0))
%                 u_nPlus1 (i,j) = max(max(U0));
%             end
            Mass_u(k+1) = Mass_u(k+1) + u_nPlus1 (i,j)*dx*dy ;
        end
    end


    %
    u_n  = u_nPlus1 ;
end
 
 Mass_err_u= Mass_u(2:nt+1) - Mass_u(1:nt )- Mass_S_half ;

  %理论上应该等于第nt层u质量加总的S的质量
% Mass = Mass_u
% Mass = 0;
%-----------limiter------------------------------------------
%    for j=1:ny
%        u_n(:,j) = limiter_2sides_final(nx,u_n(:,j),min(min(U0)),max(max(U0)));
%    end
%     for i=1:nx
%       u_n(i,:) = limiter_2sides_final(ny,u_n(i,:)',min(min(U0)),max(max(U0)))';
%     end
%---------------------------------------------------------------



end


