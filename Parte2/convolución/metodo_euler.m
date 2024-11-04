function h=metodo_euler(p,t0,y0,tf,N,dt,dim)
    clc;
    gt=zeros(dim,1);
    xg=zeros(dim,1);
    tg=zeros(dim,1);
    kg=0;
    gt(1)=1;
    tg(1)=t0;
    xg(1)=y0;

    #solución edo
    for k=1:N
      kg=dt*(gt(k)-p*xg(k));
      xg(k+1)=xg(k)+kg;
      tg(k+1)=tg(k)+dt;
      gt(k+1)=1;#función escalón
    end

    #derivada de orden 2
    h=0;
    for k=1:dim
      if(k==1)
        h(k)=(-3/(2*dt))*xg(k)+(4/(2*dt))*xg(k+1)+(-1/(2*dt))*xg(k+2);
      elseif k==dim
        h(k)=(3/(2*dt))*xg(k)+(-4/(2*dt))*xg(k-1)+(1/(2*dt))*xg(k-2);
      else
        h(k)=(-1/(dt*2))*xg(k-1)+0*xg(k)+(1/(dt*2))*xg(k+1);
      end
    end

endfunction


