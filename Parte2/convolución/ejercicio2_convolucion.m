function ejercicio2_convolucion
  clc
  p=0.5;  A=1.5;  N=200;  t0=0; tf=2*5/p; dt=(tf-t0)/N; y0=0; w=1;
  g(1)=0;
  tg=zeros(N+1);
  hh(1)=0.06;
  h=runge_kuta(p,t0,y0,tf,N,dt,N+1);
  for k=1:N
    tg(k+1)=tg(k)+dt;
    g(k+1)=A*sin(w*tg(k+1));
    hh(k+1)=exp(-p*tg(k+1));
  end
  figure(1)
  plot(tg, h,'ob', tg, hh, 'or', tg, g, 'og')
  grid on
  title("función h exacta, rk y g")

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
