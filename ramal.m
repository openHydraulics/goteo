 %Cálculo de la distribución de presión (y caudal) en un ramal de riego por goteo
function [h,q]=ramal(h0,n,k,x,CVm,le,s,D,I0,varManuf)

  vectorUnos=ones(n,1);
  matrizAcum=triu(ones(n),0);

  xR=transpose(matrizAcum)*(s.*vectorUnos);
  zR=-I0.*xR;

  h=h0.*vectorUnos;
  hant=0.*vectorUnos;

  while sum(abs(h-hant))>1e-3;
    hant=h;
    q=(k.*h.^x).*varManuf;
    h=h0.*vectorUnos-zR-transpose(matrizAcum)*(0.465.*(matrizAcum*q).^1.75.*D.^-4.75.*(le+s));
  end

end