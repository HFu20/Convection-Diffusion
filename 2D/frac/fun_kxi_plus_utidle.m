function [ kxi_plus,  Gk,count ] = fun_kxi_plus_utidle( dt, u_tidle, u_n, min_u0, max_u0 , kxi, dx,dy )

kxi_minus = 0; % 赋初值
kxi_k = 1e-8 ; %赋初值
% eps = 1e-10 ;
w = dx*dy ;  %权重
% eps =  1.0e-16; %容许误差
count = 0; %迭代次数计数器


%----------------------------------
% 先不看时间层的变化，都是在第 n+1 层迭代，
% 以通过第 n 层的 kxi 计算出新的kxi_plus (n+1层)
%-------------------------------------------
yinta_minus = dt*(   kxi_minus - kxi   ) ; % alp_k = 1;
v_minus = u_tidle +  yinta_minus ; %  原文多了dt
Gminus_1 = 0 ; %表示( k-1 )层的四个求和项（只需计算一遍）
Gminus_2 = 0 ;
Gminus_3 = 0 ;
Gminus_4 = 0;

for i = 1: length( u_tidle( :,1 ) )
    for j = 1: length( u_tidle( 1,: ) )

        if  v_minus( i, j ) <= min_u0
            Gminus_3 = Gminus_3+ min_u0* w ;
        elseif  v_minus( i, j ) >= max_u0
            Gminus_2 = Gminus_2 +max_u0* w;
        elseif v_minus( i, j )<max_u0 &&  v_minus(i,j )>min_u0
            Gminus_1 = Gminus_1 + v_minus(i,j )* w ;
        end

        Gminus_4 = Gminus_4 + u_n(i,j )* w ;

    end
end
Gminus = Gminus_1 + Gminus_2 + Gminus_3 - Gminus_4;

while abs( kxi_k - kxi_minus ) >1e-10 && count<20   %eps =  1.0e-6;

    Gk_1 = 0 ;   %表示第k层的四个求和项，需复制后迭代计算
    Gk_2 = 0 ;
    Gk_3 = 0 ;
    Gk_4 = 0;

    count = count + 1;
    %     if count > 5
    %         disp('error:不收敛');
    %         return;
    %     end

    % 需要用到 Fn( kxi_k ) = Gn( yinta_k )和 Fn( kxi_minus ) = Gn( yinta_minus )，
    % 即 Gk 和 Gminus
    yinta_k = dt*(   kxi_k - kxi  ) ; %对于时间二阶格式，B_{k-1} = 1 ;
    v_k = u_tidle + yinta_k ;
    % yinta_k 为判断所在哪个区域的自变量，v_k为判断的整体式子及求和的值
    % 给出 Gk 和 Gminus 的表达式
    for i = 1: length( u_tidle( :,1 ) )
        for j = 1: length( u_tidle( 1,: ) )
            if  v_k(i,j) <= min_u0
                Gk_3 = Gk_3+ min_u0*w ; % w =1;
            elseif  v_k(i,j ) >= max_u0
                Gk_2 = Gk_2 + max_u0* w ;
            elseif v_k(i,j )<max_u0 &&   v_k(i,j )>min_u0
                Gk_1 = Gk_1 + v_k(i,j )* w ;
            end
            Gk_4 = Gk_4 + u_n( i,j )* w ;  % （u^n,1）的内积
        end
    end
    Gk = Gk_1 + Gk_2 + Gk_3 - Gk_4;
    %     迭代求解 kxi


    %     if Gk - Gminus  ~= 0
    kxi_plus = kxi_k - Gk*( kxi_k - kxi_minus ) / ( Gk - Gminus  )  ;  % 分母=0
    kxi_minus = kxi_k ;
    kxi_k = kxi_plus;
    Gminus =  Gk;
end

end

    %% method 1
    %          if abs( Gk - Gminus )< eps
    %               kxi_plus = kxi + eps ;  % 会让kxi很大
    %          else
    %          kxi_plus = kxi_k +10*eps*sign( kxi_k - kxi_minus );
    %% method 2
    %      else   kxi_plus = kxi ;  % 会让kxi很大
    %% method 3
    %    kxi_plus = kxi_k - Gk*( kxi_k - kxi_minus ) / ( Gk - Gminus +
    %   eps )  ;  % 会让kxi很大
    %% method 4
    %     elseif  abs( min(Gk,Gminus) ) < eps*1e-2
    %          kxi_plus = kxi ;
    %     elseif Gk - Gminus  == 0
    %          kxi_plus = (kxi_minus + kxi_k)/2 ;  %+10*eps*sign( kxi_k - kxi_minus );
    %% method 5 取倒数
    %     else
    %          inv_kxi_plus = ( Gk - Gminus ) / ( Gk*kxi_minus - kxi_k*Gminus ) ;
    %%
    %         fprintf( 'Gk = Gminus= %f  at %d *dt \n ',Gk, count )
    %          return ;
    %     end




    %     end


