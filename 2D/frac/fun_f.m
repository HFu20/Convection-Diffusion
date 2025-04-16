function f = fun_f(u)

  f =  u.*u./(u.*u+(1-u).*(1-u)).*(1-5*(1-u).*(1-u));

end

