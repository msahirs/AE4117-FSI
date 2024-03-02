function [As Asf Afs Af] = matrices(N,m,k)

As  = zeros(2);
Asf = zeros(2,2*N);
Afs = zeros(2*N,2);
Af  = zeros(2*N);

As(1,2)      = -k/m;
As(2,1)      = 1;

Asf(1,N)     = 1/m;

Afs(N,1)     = -N;

Af(1:N-1    ,N+2:2*N  ) = Af(1:N-1    ,N+2:2*N  ) - eye(N-1);
Af(2:N      ,N+1:2*N-1) = Af(2:N      ,N+1:2*N-1) + eye(N-1);
Af(N+1:2*N-1,2:N      ) = Af(N+1:2*N-1,2:N      ) - eye(N-1);
Af(N+2:2*N  ,1:N-1    ) = Af(N+2:2*N  ,1:N-1    ) + eye(N-1);

Af(1,N+1) = Af(1,N+1) - 1;
Af(N+1,1) = Af(N+1,1) + 1;

Af(N,2*N) = Af(N,2*N) + 1;
Af(2*N,N) = Af(2*N,N) - 1;

Af = 0.5 * N * Af;

end