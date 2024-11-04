function ejercicio3_convolucion
  clc
  q=0.5;  A=2;  sigma=0.1;  omega=1.5;  t0=0; tf=1*max(5/sigma, 5/q);
  N=500;  dt=(tf-t0)/N; y0=0;

  tg=zeros(N+1);  hh(1)=0.06;
  h=runge_kuta_sistemas(sigma,omega,t0,y0,tf,N,dt,N+1);
  g(1)=A*exp(-q*tg(1));
  for k=1:N
    tg(k+1)=tg(k)+dt;
    g(k+1)=A*exp(-q*tg(k+1));
    hh(k+1)=(1/omega)*exp(-sigma*tg(k+1))*sin(omega*tg(k+1));
  end
  figure(1)
  plot(tg, h,'ob', tg, hh, 'or', tg, g, 'og')
  grid on

  #cálculo convolución con función de octave
  y_conv=dt*conv(hh,g);
  y_conv_rk=dt*conv(h, g);

  #cálculo convolución
  for k=1:N+1
    y_conv_for(k)=0;
    for j=1:k
      y_conv_for(k)=y_conv_for(k)+hh(k+1-j)*g(j)*dt;
    endfor
  endfor

  #solución en el dominio complejo
  h_tdf=fft(hh);
  h_rk_tdf=fft(h);
  g_tdf=fft(g);

  yc_tdf=h_tdf.*g_tdf;
  yc_rk_tdf=h_rk_tdf.*g_tdf;
  yc_inv_tdf= dt*ifft(yc_tdf);
  yc_rk_inv_tdf=dt*ifft(yc_rk_tdf);

  figure(2)
  plot(tg, y_conv(1:N+1), '-r', tg, yc_inv_tdf, '-b', tg, yc_rk_inv_tdf, '-g', tg, y_conv_rk(1:N+1), '-y');
  grid on

  figure(3)
  subplot(3,1,1)
    stem(abs(h_tdf(1:N/10)))
    grid on
    title("Módulo de h_tdf")
  subplot(3,1,2)
    stem(abs(g_tdf(1:N/10)))
    grid on
    title("Módulo de g_tdf")
  subplot(3,1,3)
    stem(abs(yc_tdf(1:N/10)))
    grid on
    title("Módulo de yc_tdf")
endfunction

function h=runge_kuta_sistemas(sigma,omega,t0,y0,tf,N,dt,dim)
    t=zeros(1,dim);
    y=zeros(2,dim);
    w=1/2;

    t(1,1)=t0; y(1,1)=y0; y(2,1)=y0;

    for k=1:N
      k1=dt*funcion(y(:,k),sigma,omega);
      tg=t(k)+dt*(1/(2*w));
      yg=y(:,k)+k1/(2*w);
      k2=dt*funcion(yg,sigma,omega);
      y(:,k+1)=(y(:,k)+(1-w)*k1+w*k2);
      t(1,k+1)=(t(1,k)+dt);
    endfor

    #derivada de orden 2
    h=0;
    for k=1:dim
      if(k==1)
        h(k)=(-3/(2*dt))*y(1,k)+(4/(2*dt))*y(1,k+1)+(-1/(2*dt))*y(1,k+2);
      elseif k==N+1
        h(k)=(3/(2*dt))*y(1,k)+(-4/(2*dt))*y(1,k-1)+(1/(2*dt))*y(1,k-2);
      else
        h(k)=(-1/(dt*2))*y(1,k-1)+0*y(1,k)+(1/(dt*2))*y(1,k+1);
      end
    endfor

endfunction

function [fy]=funcion(y,sigma, omega)
  fy(1,1)=0*y(1)+1*y(2)+0;
  fy(2,1)=-(sigma^2+omega^2)*y(1)+(-2*sigma)*y(2)+1;
endfunction
