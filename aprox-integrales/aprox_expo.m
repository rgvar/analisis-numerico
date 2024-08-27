# aprox integral función exponencial método trapecios multiple

function aprox_expo(N, p, w)
    tf = 1;
    t0 = 0;
    dt = (tf - t0) / N;


    for k=1:N+2
        tg(k) = (k-1) * dt;
        yg(k) = e^(p*tg(k)) * sin(w*tg(k));

    endfor


    It = 0;
    for k=1 : N + 1
        It = It + (yg(k) + yg(k + 1)) * dt/2;
        # fprintf("k: %d; tg=%.2f, yg=%.4f, It=%.10f\n ", k,tg(k),yg(k),It);
    endfor



    figure (1)
    plot (tg, yg, 'ob')
    grid on


    format long
    disp(It);


endfunction
