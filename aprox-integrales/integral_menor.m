function retval = integral_menor (N, tf)
    t0 = 0;
    dt = ( tf - t0 ) / N;
    Im = 0;
    for o = 1:N+1
        t(o) =(o-1) * dt;
        y(o) = sin(pi*t(o));
        if (o == 1)
            Im = 0;
        else
            Im = Im + y(o - 1);
        endif
    endfor
    Im = Im * dt;
    figure (1)
    plot (t, y, 'ob')
    grid on
    retval = Im;
endfunction
