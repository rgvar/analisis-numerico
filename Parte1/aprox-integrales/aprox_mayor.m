function retval = integral_mayor (N, tf)
    t0 = 0;
    dt = ( tf - t0 ) / N;
    IM = 0;
    for o = 1:N+1
        t(o) =(o-1) * dt;
        y(o) = sin(pi*t(o));
        if (o == 1)
            IM = 0;
        else
            IM = IM + y(o);
        endif
    endfor
    IM = IM * dt;
    figure (1)
    plot (t, y, 'ob')
    grid on
    format long
    disp((IM);
endfunction
