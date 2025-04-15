function m = minmod2( Z, epsilon, dx )

if abs(Z(1)) <= epsilon*dx^2
    m = Z(1);
else
    m = minmod(Z);
end

end


function m = minmod(Z)

s = sum(sign(Z))/length(Z);

if abs(s) == 1
    m = s*min(abs(Z));
else
    m = 0;
end

end