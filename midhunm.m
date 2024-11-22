clf;
clear all;
clear workspace;
clc;
f=9;
lamb=30/f;
rt=20000;
gama=110;
D=30;
g0=1;
G0=1;
s=1;
gamar=gama*pi/180;
delgo=0;
hpbw=65*lamb/D;
hpbwr=hpbw*(pi/180);
teta0=hpbw*0.283;
tetar0=teta0*pi/180;
tetat=0;
tmin=-hpbw;
tmax=hpbw;
tinc=.001;
teta=[tmin:tinc:tmax];
tetar=teta*pi/180;
steps=((tmax-tmin)/tinc)+1;
dh=4;
for k=1:steps
    v1(k)=(s*G0*exp((-2.776)*((teta(k)-teta0)/hpbw).^2));
    v2(k)=(s*(G0+delgo)*exp((-2.776)*((teta(k)+teta0)/hpbw).^2));
    vsum(k)=v1(k)+v2(k);
    vdiff(k)=v1(k)-v2(k);
  %  dbvsum(k)=20*log10(v1(k)+v2(k));
   % dbvdiff(k)=20*log10(v1(k)-v2(k));
    ve(k)=real(vdiff(k)/vsum(k));
end
ve1=round(ve*10);
ve2=ve1/10;
lmin=50;
linc=25;
lmax=100;
lsteps=(lmax-lmin)/linc+1;
l=[lmin:linc:lmax];
for l1=1:lsteps
    if gamar==90*pi/180
        gamar=89.98*pi/180;
    elseif gamar==140*pi/180
       gamar=139.98*pi/180;
    elseif gamar==30*pi/180
       gamar=29.98*pi/180;
   else gamar=gamar;
   end
    rd(l1)=sqrt(rt*rt+l(l1)*l(l1)-2*rt*l(l1)*cos(pi-gamar));
    temp(l1)=(l(l1)/rd(l1))*sin(pi-gamar);
    tetadr(l1)=asin(temp(l1));
    tetad(l1)=tetadr(l1)*180/pi;
end
s=1;
Jmin=0;
Jmax=5;
Jinc=0.1;
J=[Jmin:Jinc:Jmax];
Jsteps=(Jmax-Jmin)/Jinc+1;
n=Jsteps;
randn('seed');
na=randn(1,n);
snrdb=0;
snr=10^(snrdb/20);
for l1=1:lsteps
    for m=1:Jsteps
        path(l1)=rt-rd(l1);
        psi(l1)=2*pi*path(l1)/(lamb/100);
        v1s(m,l1)=sqrt(s*g0*exp((-2.776)*((tetat-teta0)/hpbw)^2))+na(l1)/snr;
        %v1s(m,l1)=(s*g0*exp((-2.776)*((tetat-teta0)/hpbw)^2));
        v1j(m,l1)=(J(m)*g0*exp((-2.776)*((tetad(l1)-teta0)/hpbw)^2))*exp(j*psi(l1));
        v1(m,l1)=v1s(m,l1)+v1j(m,l1);
        v2S(m,l1)=(s*g0*exp((-2.776)*((tetat+teta0)/hpbw)^2));
        v2j(m,l1)=(J(m)*g0*exp((-2.776)*((tetad(l1)+teta0)/hpbw)^2))*exp(j*psi(l1));
        v2(m,l1)=v2S(m,l1)+v2j(m,l1);
        vsum(m,l1)=v1(m,l1)+v2(m,l1);
        vdiff(m,l1)=v1(m,l1)-v2(m,l1);
        volerr(m,l1)=real(vdiff(m,l1)/vsum(m,l1));
    end
end
volerr1=round(volerr*10);
volerr2=volerr1/10;
for l1=1:lsteps
    for m=1:Jsteps
        for k=1:steps
            if(abs(volerr(m,l1)-ve(k)))<=0.001
                angerr(m)=teta(k);
            else
            end
        end
 md(m)=abs(rt*angerr(m)*(pi/180));
    end
    hold on;
    figure(1);
    plot(J,md);
    xlabel('j/s ratio.............>');
    ylabel('miss distance in  meters...>');
    title('gamma=110');
    if l1==1
        gtext('l=50');
    elseif l1==2
        gtext('l=75');
    elseif l1==3
        gtext('l=100');
%    elseif l1==4
 %       gtext('l=500');
    else
    end
    hold off;
    grid on;
end