function h=runge_kuta(p,t0,y0,tf,N,dt,dim)
    gt=zeros(dim,1);
    y=zeros(dim,1);
    t=zeros(dim,1);
    kg=0;
    gt(1)=1;
    t(1)=t0;
    y(1)=y0;
    w=1/2;

    #solución edo
    #for k=1:N
     # kg=dt*(gt(k)-p*xg(k));
      #xg(k+1)=xg(k)+kg;
      #tg(k+1)=tg(k)+dt;
      #gt(k+1)=1;#función escalón
    #end

    for k=1:N
      k1=dt*((gt(k)-p*y(k)));
      tg=t(k)+(dt/(2*w));
      yg=y(k)+(k1/(2*w));
      k2=dt*((gt(k)-p*yg));
      y(k+1)=y(k)+(1-w)*k1+w*k2;
      t(k+1)=t(k)+dt;
      gt(k+1)=1;
    endfor

    #derivada de orden 2
    h=0;
    for k=1:dim
      if(k==1)
        h(k)=(-3/(2*dt))*y(k,:)+(4/(2*dt))*y(k+1,:)+(-1/(2*dt))*y(k+2,:);
      elseif k==dim
        h(k)=(3/(2*dt))*y(k,:)+(-4/(2*dt))*y(k-1,:)+(1/(2*dt))*y(k-2,:);
      else
        h(k)=(-1/(dt*2))*y(k-1,:)+0*y(k,:)+(1/(dt*2))*y(k+1,:);
      end
    end

endfunction
