# aproximación integral por método simpson
# función y= sin(pi*t(o))

function aprox_simpson(N)
    if mod(N,2) != 0
        return;
    endif
    t0 = 0;
    tf = 0.5;
    dt = (tf - t0) / N;
    for k=1: N + 1
        tg(k) = (k-1) * dt;
        yg(k) = sin(pi*tg(k));
    endfor
    ISim = 0;
    for k=2: 2:N
        ISim = ISim + dt*(yg(k-1)+4*yg(k)+yg(k+1))/3;
    endfor
    format long
    disp(ISim);
endfunction
