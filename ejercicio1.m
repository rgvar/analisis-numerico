function retval = ejercicio1 (N, tf)
    t0 = 0;
    dt = (tf - t0) / N;
    It = 0;
    for k=1:N+1
        t(k) = (k-1) * dt;
        y(k) = sin(t(k));
        if (k == 1)
            It = 0;
        else
            It = It + (y(k-1) + y(k)) / 2;
        endif
    endfor
    It = It*dt;
    retval = It;

endfunction
