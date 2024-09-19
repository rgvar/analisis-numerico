

function df_dospuntos_seno(x = pi/10, h = 0.1)

    dfexacta = cos(x);

    for k=1:5
    printf("h = %.5f \n", h);
    df = (sin(x + h) - sin(x)) / h;
    format long
    df
    Er = dfexacta - df;
    format long
    Er
    h = h/10;
    disp(" ");
    endfor

endfunction
