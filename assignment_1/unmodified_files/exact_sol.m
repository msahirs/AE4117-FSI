function [W] = exact_sol(o,N,t)
  x              = 0.5/N:1/N:(1-0.5/N);
  W              = zeros(2*N+2,length(t));
  W(1:2,:)       = [ -o*sin(o*t)
                        cos(o*t) ];
  W(3:N+2,:)     = - o / sin(o) * cos(o*x') * cos(o*t);
  W(N+3:2*N+2,:) = - o / sin(o) * sin(o*x') * sin(o*t);
end
