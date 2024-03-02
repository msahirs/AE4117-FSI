function [W_f W_s] = exact_sol(o,N,t)
  W_s          = [ -o*sin(o*t) ; cos(o*t) ];
  W_f          = zeros(2*N,1);
  x            = 0.5/N:1/N:(1-0.5/N);
  W_f(1:N)     = - o / sin(o) * cos(o*x') * cos(o*t);
  W_f(N+1:2*N) = - o / sin(o) * sin(o*x') * sin(o*t);
end
