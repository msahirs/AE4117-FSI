function omega = exact_omega(m,k)
  nmax_iter = 100;
  phase_min = 1.e-10;
  phase_max = pi;
  npi       = 0;
  phase     = 0.5 * pi;
  o_old     = -1;
  o         =  1;


  i = 0;
  while ( ~ ((o == o_old) || (i > nmax_iter)) )
      o_old = o;
      i = i + 1;
      o = npi * pi + phase_min;
      err_min = ((m * o*o - k) * sin(phase_min) - o * cos(phase_min)) / ...
                (k + o + m * o*o);
      o = npi * pi + phase_max;
      err_max = ((m * o*o - k) * sin(phase_max) - o * cos(phase_max)) / ...
                (k + o + m * o*o);
      o = npi * pi + phase;
      err = ((m * o*o - k) * sin(phase) - o * cos(phase)) / ...
            (k + o + m * o*o);
      if (err_min * err_max > 0)
        baderror
      end
      if (err * err_min < 0) 
        phase_max = phase;
        phase = phase_min + (phase - phase_min) / (err_min - err) * err_min;
      else
        phase_min = phase;
        phase = phase + (phase_max - phase) / (err - err_max) * err;
      end
  end
  if ( o_old ~= o )
    disp('Could not find correct mode...');
  end
  omega = o;
end
