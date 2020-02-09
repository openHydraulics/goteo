nr=10;ng=20;h0=10;k=1.25;x=0.54;CVm=0.05;le=0.5;s=0.5;D=13.6;I0=0.05;

h=h0.*ones(ng,nr);
q=k.*h.^x;

%Caudal en los tramos del ramal
Qr=triu(ones(ng),0)*q

%Caudal en los tramos de la tubeŕia portarramal
Qpr=(ones(1,ng)*q)*tril(ones(nr),0)