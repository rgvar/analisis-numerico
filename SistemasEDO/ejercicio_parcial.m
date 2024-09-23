function ejercicio_parcial

    # y1' = 0 * y1 + 1 * y2 + 0 *sin(3*t)
    # y2' = -(2^2) * y1 - 0.02 * y2 - 1 * sin(3*t)
    # x1(0) = 0
    # x2(0) = 0

    tf = 12.5664;
    t0 = 0;
    NDt = 3000;
    Dt = (tf-t0) / NDt;

    t = zeros(1,NDt);
    y = zeros(2,NDt);
    k1 = zeros(2,1);

    f = @(y1,y2,t) [ y2 ; -4*y1 - 0.02*y2 - sin(3*t)];


    # EULER EXPLÍCITO
    t(1) = 0;
    y(:,1) = 0;
    for j=1 : NDt-1
        k1 = Dt * f( y(1,j) , y(2,j) , t(j) );
        t(j+1) = t(j) + Dt;
        y(:,j+1) = y(:,j) + k1;
    endfor

    fprintf("y1(NDt)= %.4f \ny2(NDt)= %.4f \n",y(1,NDt),y(2,NDt))

    figure(1)
    plot(t,y(1,:),'r')
    grid on

    # TRAPECIOS COMPUESTOS
    It=0;
    for j=1:NDt - 1
        It = It + ( Dt * ( y(1,j)^2 + y(1,j+1)^2 )) / 2;
    endfor

    f2= @(y1,y2) y1 * y2;

    It2=0;
    for j=1:NDt - 1
        It2=It2 + Dt * (f2(y(1,j+1),y(2,j+1)) + f2(y(1,j),y(2,j))) / 2;
    endfor

    It
    It2


endfunction
