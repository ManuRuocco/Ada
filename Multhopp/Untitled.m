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
    gamma=zeros(1,mmax);
    gamb=zeros(1,mmax);
    gama=zeros(1,mmax);
    disgam=zeros(1,mmax);
    gam=zeros(1,mmax);
    a=zeros(mmax);
    
    dummy=pi/(m+1);
    mmaxm=mmax-1;

