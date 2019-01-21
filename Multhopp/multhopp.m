
function multhopp
clc 
clear all

%INPUT%%%%%%%%%%%%%%%%%%
lab=fopen('In','r');
fgetl(lab);
ie=str2num(fgetl(lab));
fgetl(lab);
ca=str2num(fgetl(lab));
fgetl(lab);
eta=str2num(fgetl(lab));
eta1=eta(1);eta2=eta(2);
fgetl(lab);
aerdat=str2num(fgetl(lab));
alam=aerdat(1);ar=aerdat(2);er=aerdat(3);ee=aerdat(4);
fgetl(lab);
rprt=str2num(fgetl(lab));
clar=rprt(1);clae=rprt(2);
fgetl(lab);
dum=str2num(fgetl(lab));
m=dum(1);alfa=dum(2);
fclose(lab);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (ca==0)  %senza alettoni
    eta1=1.1;
    eta2=1.1;
end
mmax=(m+1)/2; 
%    CALCOLO FATTORE CORRETIVO FORMULA DI PRANDTL***********************      
if(ar>3)
    enne=.5+(4./(pi^2)/ar)*(log(pi*ar)-(7./8.));
else
    enne=1.-.25*ar+.067147*(ar^2)-.0062767*(ar^3);
end                 
%    AZZERAMENTO MATRICI E VETTORI**************************************
for i=1:mmax
    gamma(i)=0;
    gamma(m+1-i)=0;
    gamb(i)=0;
    gama(i)=0;
    disgam(i)=0;
    gam(i)=0;
    for j=1:mmax
        a(i,j)=0;
    end
end
dummy=pi/(m+1);
mmaxm=mmax-1;
%    CALCOLO SENI E COSENI**********************************************
for i=1:mmaxm
    mp1mi=m+1-i;
    teta=i*dummy;
    s1=sin(teta);
    c1=cos(teta);
    seno(i)=s1;
    seno(mp1mi)=s1;
    coseno(i)=c1;
    coseno(mp1mi)=-c1;
end
seno(mmax)=1.0;
coseno(mmax)=0.0;
%    CALCOLO VALORI DI SVERGOLAMENTO E CLA PER ALI AVENTI VARIAZIONI LINEARI
%    DI TALI PARAMETRI(o paraboliche)***************************************      
for i=1:mmax
    e(i)=er+(ee-er)*coseno(i);
    cla(i)=clar+(clae-clar)*coseno(i);
    %      e(i)=er+(ee-er)*coseno(i)
    %      e(i)=ee*coseno(i)**2
end
dummy=(m+1)/4;
dummy2=pi./180./2/enne;
conver=pi/180/2./enne;
is=0;
ib=0;
alfar=0.;
test=1;
% termine noto alfa=1 grado
while test
    if(is==1), mmax=mmax-1; end %41
    %    COSTRUZIONE DELLLA MATRICE D'INFLUENZA*****************************
    for i=1:mmax
        si=seno(i);
        ci=coseno(i);
        for j=1:mmax
            sj=seno(j);
            cj=coseno(j);
            if (i==j)
                if (ie==0), s(i)=ar*(1.+alam)/(1.-ci*(1.-alam))/cla(i); end
                if (ie==1), s(i)=pi*ar/(2.*sqrt(1.-ci*ci))/cla(i); end
                a(i,j)=((m+1)/4./si)+s(i)/2./enne;
                continue
            end
            if(mod(i-j,2)==0),  continue; end
            b1=sj/(cj-ci)^2/(m+1);
            b2=sj/(-cj-ci)^2/(m+1);
            %    SUDDIVISIONE IN CARICO SIMMETRICO E IN CARICO ANTISIMMETRICO*******
            if(is==1)
                a(i,j)=-b1+b2;
            else
                a(i,j)=-b1-b2;
                if(j==mmax),  a(i,j)=a(i,j)/2; end
            end
        end
    end
    mmaxp=mmax+1;
    mmaxq=mmaxp+1;
    dalfar=conver;
    if(ib~=1)
        %    CARICO BASICO******************************************************
        for k=1:mmax
            a(k,mmaxp)=1; 
            a(mmaxp,k)=seno(k)*2;
            a(k,mmaxq)=e(k)*conver;
        end
        a(mmaxp,mmax)=1;   
        a(mmaxp,mmaxp)=0;
        a(mmaxp,mmaxq)=0;
        m1=mmaxp;
        m2=mmaxq;
    else
        m1=mmax;
        m2=mmaxp;
        if(is~=1) 
            %    TERMINE NOTO CASO SIMMETRICO***************************************
            for i=1:mmax
                a(i,mmaxp)=dalfar;
            end
        else
            %    TERMINE NOTO CASO ANTISIMMETRICO***********************************
            for i=1:mmax
                ci=coseno(i);
                a(i,mmaxp)=0.;
                if((ci>=eta1)&(ci<=eta2)), a(i,mmaxp)=dalfar; end
            end
            %    SOLUZIONE DEL SISTEMA**********************************************
        end
    end
    test=0;a1=a;
    a=ssystem(m1,m2,cla,alam,ar,a1);
    if (ischar(a))
        uiwait(msgbox('DETERMINANTE NULLO','ERRROR!'));
        return
    end
    if(ib~=1)
        %    gamma simmetrica basica********************************************  
        for i=1:mmax
            gamb(i)=a(i,mmaxq);
        end
        alfar=a(mmaxp,mmaxq);
        alfar=alfar*180/pi;
        ib=1;
        test=1;
        %        continue
    else
        if(is~=1)
            %    gamma simmetrica addizionale***************************************
            for i=1:mmax
                gama(i)=a(i,mmaxp);
            end
            is=1;
            test=1;
            %            continue
        else
            %     gamma antisimmetrica*********************************************  
            for i=1:mmax
                disgam(i)=a(i,mmaxp);
                %        disgam(i)=0          
            end
        end
    end
end
mmax=mmaxp;
%    INCREMENTO LINEARE DELL'INCIDENZA**********************************
for i=1:(m+1)/2
    gama(i)=gama(i)*(alfa+alfar);
    disgam(i)=disgam(i)*(alfa+alfar);
end
%    GAMMA TOTALE*******************************************************
for i=1:(m+1)/2
    gam(i)=gamb(i)+gama(i);
    mp1mi=m+1-i;
    gamma(i)=gam(i)+disgam(i);
    gamma(mp1mi)=gam(i)-disgam(i);
end
%    CALCOLO DI CDI e di cl******************************************************   

dummy1=0;
dummy2=0;
% for i=1:(m+1)/2
%     ai(i)=(-gamma(i)*s(i)+(alfa+e(i))*pi/180)/enne/2.;
% end
% for i=1:(m-1)/2
%     dummy1=dummy1+gamma(i)*seno(i);
%     dummy2=gamma(i)*seno(i)*ai(i)+dummy2;
% end
% cl=2*dummy1*pi*ar/(m+1)+gamma((m+1)/2)*pi*ar/(m+1);
% cd=2*dummy2*pi*ar/(m+1)+gamma((m+1)/2)*pi*ar*ai((m+1)/2)/(m+1);
for i=1:m
    dummy1=dummy1+gamma(i)*seno(i);
    dummy2=gamma(i)*seno(i)*coseno(i)+dummy2;
end
cl=dummy1*pi*ar/(m+1);
cm=dummy2*pi*ar/(m+1)/2; %funziona!
dummy1=0;
for i=1:m
    
    s1=seno(i);
    c1=coseno(i);
    dummy2=0.0;
    for j=1:m
        if((i==j)|(mod(i-j,2)==0)),  continue; end
        sj=seno(j);
        cj=coseno(j);
        dummy2=gamma(j)*sj/(cj-c1)^2/(m+1)+dummy2;
    end
    c=dummy/s1*gamma(i);
    dummy1=dummy1+gamma(i)*s1*(c-dummy2);
    
end
cd=dummy1*pi*ar/(m+1);
%end
azl=-alfar,cl,cd,cm
%    output su file%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out1=[coseno(1:(m+1)/2);gamb(1:(m+1)/2);gama(1:(m+1)/2);gamma(1:(m+1)/2);disgam(1:(m+1)/2)]';
save multhopp1.dat out1 -ascii
out2=[coseno(1:m);gamma(1:m)]';
save multhopp2.dat out2 -ascii
figure(1),plot(coseno(1:(m+1)/2),gamb(1:(m+1)/2),'b'),ylabel('carico basico'),xlabel('y/(b/2)'),grid
figure(2),plot(coseno(1:(m+1)/2),gama(1:(m+1)/2),'b'),ylabel('carico addizionale'),xlabel('y/(b/2)'),grid
figure(3),plot(coseno(1:(m+1)/2),disgam(1:(m+1)/2),'b'),ylabel('carico antisimmetrico'),xlabel('y/(b/2)'),grid
figure(4),plot(coseno(1:m),gamma(1:m),'b'),ylabel('carico totale'),xlabel('y/b'),grid
%    FINE PROGRAMMA-----------------------------------------------------


%    SUBROUTINE*********************************************************
function a=ssystem(n,m,cla,alam,ar,a1)

for i=1:n 
    pivot=i;
    xax=abs(a1(i,i));
    for k=i:n 
        if (abs(a1(k,i))-xax)>0
            pivot=k;
            xax=abs(a1(k,i));
        end
    end
    if (xax) <=0 
        %---------------->determinante nullo
        a='out';
        return
    end
    k=pivot;    
    for j=1:m
        pivot=a1(i,j);
        a1(i,j)=a1(k,j);
        a1(k,j)=pivot;
    end
    pivot=a1(i,i);
    for k=1:m
        a1(i,k)=a1(i,k)/pivot;
    end
    for k=1:n
        pivot=a1(k,i);
        if (k-i)
            for j=1:m
                a1(k,j)=a1(k,j)-pivot*a1(i,j);
            end
        end
    end
end
a=a1;