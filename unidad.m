nR=25;ng=100;h0=10;k=1.25;x=0.54;CVm=0.05;leg=0.75;sr=0.5;Dr=13.6;I0r=0.0;lepR=1.5;spR=1;DpR=50;I0pR=0;

%Se calcula la distribución de presión hpR en la tubería portarramal
%se utiliza el paquete Image -> pkg load image
pkg load image;
tic;

vectorUnos=ones(1,nR);
matrizAcum=triu(ones(nR),0);

xpR=(spR.*vectorUnos)*matrizAcum;
zpR=-I0pR.*xpR;

hpR=h0.*vectorUnos;
hpRant=0.*vectorUnos;

%Variación manufactura
varManuf=1+CVm.*randn(ng,nR);

cont=0;

while max(abs(hpR-hpRant))>1e-3;
  hpRant=hpR;
  for j=nR:-1:1
    %Presiones y caudales en los goteros del ramal j
    [h,q]=ramal(hpRant(j),ng,k,x,CVm,leg,sr,Dr,I0r,varManuf(:,j));
    hr(:,nR+1-j)=h;%Matriz de presiones en los goteros de la unidad
    qr(:,nR+1-j)=q;%Matriz de caudales de los goteros de la unidad
  end

  %Caudal en los tramos de la tubeŕia portarramal
  QpR=(ones(1,ng)*qr)*transpose(matrizAcum);
  hpR=h0.*vectorUnos-zpR-(0.465.*QpR.^1.75.*DpR.^-4.75.*(lepR+spR))*matrizAcum;
  ++cont;
end

toc

qh=(k.*hr.^x);

%Resultados
CVq=std2(qr)/mean2(qr)
CVqh=std2(qh)/mean2(qh)
CVqm=sqrt(CVq^2-CVqh^2)

surf(hr)
hold on
surf(qr)