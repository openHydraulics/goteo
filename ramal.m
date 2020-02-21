 %Cálculo de la distribución de presión (y caudal) en un ramal de riego por goteo
 %Devuelve la distribución de presión y caudales de los emisores
function [h,q]=ramal(h0,n,k,x,CVm,le,s,D,I0,varManuf)

  vectorUnos=ones(n,1);
  matrizAcum=triu(ones(n),0);

  %Distancia al origen y cota de cada emisor
  xR=transpose(matrizAcum)*(s.*vectorUnos);
  zR=-I0.*xR;

  h=h0.*vectorUnos;
  hant=0.*vectorUnos;

  %Cálculo iterativo de la distribución de presiones en el ramal.
  while max(abs(h-hant))>1e-3;
    hant=h;
    q=(k.*h.^x).*varManuf;
    q=q.*(q>0);
    %Las pérdidas de carga se calculan con Blasius, viscosidad a 20ºC y longitud equivalente en las inserciones de ramal
    h=h0.*vectorUnos-zR-transpose(matrizAcum)*(0.465.*(matrizAcum*q).^1.75.*D.^-4.75.*(le+s));
  end

end