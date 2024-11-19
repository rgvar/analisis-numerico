function GLOBAL

    g = load('registro_231121.txt');
    N = length(g);
    s = 1.17;
    dt = 0.02;
    tf = (N-1)*dt;

    t(1) = 0;

    for i=1 : N
        t(i) = (i-1)*dt;
        h(i) = exp(-s*t(i));
    end

    h2 = dt*conv(h,h);
    h3 = dt*conv(h,h2);

    yc = dt*conv(h3,g);

    figure(1);
    title("g");
    plot(t,g,'b');
    grid on;

    #Función escalón unitaria
    for i=1:N
        gu(i)=1;
    end

    y(1,1) = 0;
    y(2,1) = 0;
    y(3,1) = 0;

    #EDO
    for i=1: N-1
        k1 = dt * f_pend(g(i), y(:,i));
        y(:,i+1) = y(:,i) + k1;
    end

    figure(2)
    plot(t, y(3,:) ,'b');
    grid on;

    figure(3);
    plot(t, h, 'b');
    grid on;

    figure(4);
    plot(t, h3(1:N),'b');

    figure(5);
    plot(t, yc(1:N), 'b');
    grid on;

    #Trapecios
    X32 = y(3,:);
    Xerr = abs(y(3,:) - yc(1:N));

    for i=1 : N
        X32_cuad(i) = X32(i)^2;
        Xerr_cuad(i) = Xerr(i)^2;
    endfor

    #I1 = 0;
    #for i=1 : N-1
    #    I1 = I1 + (dt * ((y(3,i) + y(3,i+1))/2)^2);
    #end
    #
    #I1
    #
    #I2 = 0;
    #for i=1 : N-1
    #    I2 = I2 + (dt * (((y(3,i) + y(3,i+1))/2) - ((yc(i)+yc(i+1))/2))^2);
    #end
    #
    #Id = I2

    I1 = 0;
    for i=1 : N-1
        I1 = I1 + (dt * ((X32_cuad(i) + X32_cuad(i+1))/2));
    end
    I1

    I2 = 0;
    for i=1 : N-1
        I2 = I2 + (dt * ((Xerr_cuad(i) + Xerr_cuad(i+1))/2));
    end

    Id = I2

    figure(6);
    plot(t, yc(1:N), 'b', t, y(3,:), 'r');
    grid on;

end

function [fy] = f_pend(x,z)

    s = 1.17;

    fy(1,1) = -s * z(1) + x(1);
    fy(2,1) = 1*z(1) + (-s*z(2)) + 0;
    fy(3,1) = 1*z(2) + (-s*z(3)) + 0;


end
