function ejercicio1_convolucion
  #consultar por valores iniciales entre la exacta y la obtenida con rk
  #diferencias entre la obtenida por convolución y por tdf
  clc
  p=0.5;  A=0.2;  N=200;  t0=0; tf=2*5/p; dt=(tf-t0)/N; y0=0; g(1)=A*1;
  tg=zeros(N+1);
  hh(1)=0.06;
  h_euler=metodo_euler(p,t0,y0,tf,N,dt,N+1);
  h_rk=runge_kuta(p,t0,y0,tf,N,dt,N+1);
  for k=1:N
    tg(k+1)=tg(k)+dt;
    g(k+1)=A*1;
    hh(k+1)=exp(-p*tg(k+1));
  end
  #h_euler(1:5)
  #h_rk(1:5)
  #hh(1:5)
  #tg(1:5)

  figure(1)
  plot(tg, h_euler,'ob', tg, hh, 'xr', tg, h_rk, 'og')
  grid on
  title("función h obtenida con rk, euler adelante y exacta")

  figure(2)
  plot(tg, hh, 'or', tg, g, 'ob')
  grid on
  title("función h exacta y g")

  #cálculo convolución con función de octave
  y_conv_hh=dt*conv(hh,g);
  y_conv_hrk=dt*conv(h_rk,g);

  #cálculo convolución
  for k=1:N+1
    y_conv_for(k)=0;
    for j=1:k
      y_conv_for(k)=y_conv_for(k)+hh(k+1-j)*g(j)*dt;
    endfor
  endfor

  y_conv_exacta(1)=0;
  for k=1:N
    y_conv_exacta(k+1)=(A/p)*(1-exp(-p*tg(k+1)));
  endfor

  figure(3)
  plot(tg, y_conv_hh(1:N+1), '-r', tg, y_conv_exacta, '-b', tg, y_conv_for, '-g', tg, y_conv_hrk(1:N+1), '-y');
  grid on
  title("Resultados convolución a partir de h exacta(roja), con for(verde), exacta(azul) y rk(amarilla).\n La roja coincide con la verde");

  figure(4)
  plot(tg, y_conv_hh(1:N+1), '-r', tg, y_conv_hrk(1:N+1), '-b');
  grid on
  title("Resultados convolución exacta(roja) y rk(azul)")

  #convolución con TDF
  hh_tdf=fft(hh);
  hrk_tdf=fft(h_rk);
  g_tdf=fft(g);

  #hh_tdf(1:5)
  #hrk_tdf(1:5)

  yc_tdf=hh_tdf.*g_tdf;
  yc_rk_tdf=hrk_tdf.*g_tdf;

  for k=1:N+1
    y_tdf(k)=hh_tdf(k)*g_tdf(k);
    y_rk_tdf(k)=hrk_tdf(k)*g_tdf(k);
  endfor

  yc_inv_tdf=dt*ifft(yc_tdf);
  yc_rk_inv_tdf=dt*ifft(yc_rk_tdf);

  y_inv_tdf=dt*ifft(y_tdf);
  y_rk_inv_tdf=dt*ifft(y_rk_tdf);

  #yc_inv_tdf(1:5)
  #yc_rk_inv_tdf(1:5)
  #y_conv_hh(1:5)
  #y_conv_hrk(1:5)

  figure(5)
  plot(tg, y_conv_hh(1:N+1), '-r', tg, y_conv_hrk(1:N+1), '-b', tg, yc_inv_tdf, '-g', tg, yc_rk_inv_tdf, '-y');
  grid on
  title("Resultados convolución exacta(roja) y rk(azul), tdf con h exacta(verde) y tdf con rk(amarilla)")

  figure(6)
  plot(tg, y_conv_hh(1:N+1), '-r', tg, y_conv_hrk(1:N+1), '-b', tg, y_inv_tdf, '-g', tg, y_rk_inv_tdf, '-y');
  grid on
  title("Resultados convolución exacta(roja) y rk(azul), tdf con h exacta(verde) y tdf con rk(amarilla)")

endfunction


