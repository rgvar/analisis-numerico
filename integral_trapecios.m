function retval = integral_trapecios (N, tf)
    t0 = 0;
    dt = ( tf - t0 ) / N;
    It = 0;
    for o = 1:N+1
        t(o) =(o-1) * dt;
        y(o) = sin(pi*t(o));
        if (o == 1)
            It = 0;
        else
            It = It + (y(o-1)+ y(o))/2;
        endif
    endfor
    It = It * dt;
    figure (1)
    plot (t, y, 'ob')
    grid on
    retval = It;
endfunction
