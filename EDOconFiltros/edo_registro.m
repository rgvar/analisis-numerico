function edo_registro
    fs=16000;
    Dt=1/fs;
    p=12;

    arreglold= load("Registro_240901.txt", "-ascii");
    NDt = length(arreglold);

    t = zeros(NDt, 1);
    y = zeros(NDt, 1);
    h = zeros(NDt, 1);

    t(1) = 0;
    y(1) = 0;
    h(1) = arreglold(1);
    for j=1: NDt-1
        k1=Dt * (arreglold(j) - p*y(j));
        t(j+1) = t(j) + Dt;
        y(j+1) = y(j) + k1;
        h(j+1) = arreglold(j);
    endfor

    figure(1)
    plot(t,y,'r')
    grid on

endfunction
