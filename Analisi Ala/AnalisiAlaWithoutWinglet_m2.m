clc; close all; clear all;
format long

b_wo_wl=58.3;
root_kinkAF=load('sc(2)-0614chiuso.dat');
tipkAF=load('sc(2)-0610chiuso.dat');

c_root_0=11.982; xLE_root_0=18.626; yLE_root_0=0.0;   zLE_root_0=-1.464;
c_root_f=10.500; xLE_root_f=20.400; yLE_root_f=2.820; zLE_root_f=-1.15;
c_kink=7.200;      xLE_kink=24.35;    yLE_kink=9.1;     zLE_kink=-0.45;
c_tip=2.57;        xLE_tip=36.973;    yLE_tip=29.15;     zLE_tip=0.86;
c_wl_tip=0.93;   xLE_wl_tip=39.427; yLE_wl_tip=30.15;  zLE_wl_tip=2.161;

aol_Root =-3.8249*(pi/180);
aol_tip = -3.5720*(pi/180);

Cl_slope_root=6.64    
Cl_slope_tip=6.56

twist_root = 0*(pi/180);
twist_tip = -4*(pi/180);

sweep_le = atan((xLE_tip-xLE_root_0)/(b_wo_wl/2)) *(180)/pi
A1=(c_root_0+c_kink)*yLE_kink*0.5;            %Aerea primo pannello
A2=(c_kink+c_tip)*(yLE_tip-yLE_kink)*0.5;     %Aerea primo pannello
S_recalc=2*(A1+A2)
c0=4*S_recalc/(pi*b_wo_wl);

y=linspace(0,b_wo_wl/2,1000);
c_ell=c0*(2/b_wo_wl)*sqrt(0.25*b_wo_wl^2-y.^2);

S_ell=2*trapz(y,c_ell)

c_vs_y=nan(1,length(y));
c_equiv_vs_y=nan(1,length(y));
Cl_slope_vs_y=nan(1,length(y));
aol_vs_y=nan(1,length(y));
twist_vs_y=nan(1,length(y));
airfoil_1=nan(size(root_kinkAF,1),2,length(y));
airfoil_2=nan(size(tipkAF,1),2,length(y));
        
mean_Cl_slope=(2/S_recalc)*(Cl_slope_root*A1+Cl_slope_tip*A2);

c_root_equiv=S_recalc*2/b_wo_wl-c_tip
rast_equiv=c_tip/c_root_equiv
AR=b_wo_wl^2/S_recalc

for i=1:length(y)
c_vs_y(i)=c(y(i),c_root_0,c_kink,c_tip,yLE_kink,yLE_tip,1);
c_equiv_vs_y(i)=c(y(i),c_root_equiv,c_kink,c_tip,yLE_kink,yLE_tip,0);

 if y(i)<yLE_kink
    airfoil_1(:,:,i) = root_kinkAF*c_vs_y(i);
    Cl_slope_vs_y(i) = Cl_slope_root;
    aol_vs_y(i) =aol_Root;
    twist_vs_y(i) = twist_root;
 else
    airfoil_2(:,:,i)=tipkAF*c_vs_y(i);
    Cl_slope_vs_y(i)=Cl_slope_root+((Cl_slope_tip-Cl_slope_root)/(yLE_tip-yLE_kink))*(y(i)-yLE_kink);
    aol_vs_y(i)=aol_Root + ((aol_tip-aol_Root)/(yLE_tip-yLE_kink))*(y(i)-yLE_kink);
    twist_vs_y(i)=twist_root + ((twist_tip-twist_root)/(yLE_tip-yLE_kink))*(y(i)-yLE_kink);
 end
end

AOL_rad=(2/S_recalc).*trapz(y,c_vs_y.*(aol_vs_y-twist_vs_y))

k=Cl_slope_vs_y/mean_Cl_slope;

xLE_vs_y=xLE_root_f+((xLE_tip-xLE_root_f)/(yLE_tip-yLE_root_f))*(y-yLE_root_f);
xTE_vs_y=xLE_vs_y+c_vs_y;
xTE_equiv_vs_y=xLE_vs_y+c_equiv_vs_y;

mac=2*trapz(y,c_vs_y.^2)/S_recalc;

xLE_mac=2*trapz(y,c_vs_y.*xLE_vs_y)/S_recalc;

yLE_mac=2*trapz(y,c_vs_y.*y)/S_recalc;

%Calcolo del cCl addizionale per CL=1

cmean = 0.5*(c_vs_y+c_ell);
cCl_a1 = (c_ell + c_vs_y.*k)/2;
gamma_a1 = cCl_a1/(2*(b_wo_wl));
gamma_a1(end)=0;
Cl_a1=cCl_a1./c_vs_y.*k;
Cl_04= Cl_a1*0.4;
CL=trapz(y/(b_wo_wl/2),Cl_a1);

%Effetto freccia Pope ed Haney 
Cl_sweep = Cl_a1./k - ((1-(y/(b_wo_wl/2)))*(2*(1-cos(sweep_le*pi/180))));
CL_sweep=trapz(y/(b_wo_wl/2),Cl_sweep);
Cl_sweep_correct = Cl_sweep./CL_sweep.*CL;
CL_sweep_correct = trapz(y/(b_wo_wl/2),Cl_sweep_correct);
gamma_sweep = Cl_sweep.*c_vs_y/(2*b_wo_wl);
gamma_sweep_correct = Cl_sweep_correct.*c_vs_y/(2*b_wo_wl);


%carico basico

cCl_b=0.5*(c_vs_y).*Cl_slope_vs_y.*(twist_vs_y-(AOL_rad));
cCl_b(end)=0;
Cl_b=cCl_b./c_vs_y;
CL=trapz(y/(b_wo_wl/2),Cl_b)


% Cl totale
cCl_tot = cCl_b + cCl_a1;
Cl_tot = cCl_tot./c_vs_y;


subplot(1,3,1)
plot(y,c_vs_y,y,c_ell,y,0.5*(c_vs_y+c_ell),[yLE_mac yLE_mac],[0 mac],'r')
grid minor
title('Distribuzione di corde','Interpreter','latex')
xlabel('$y$ \qquad','Interpreter','latex')
ylabel('$c(y)$','Interpreter','latex')
legend('c','c ell', '\gamma_1')
set(get(gca,'ylabel'),'rotation',0)
axis equal;
axis([0 b_wo_wl/2 0 c_root_0+1])

subplot(1,3,2)
plot(y,xTE_vs_y,y,xLE_vs_y,y,xTE_equiv_vs_y,[yLE_mac yLE_mac],[xLE_mac xLE_mac+mac],'r')
axis ij
grid minor
title('Pianta semiala','Interpreter','latex')
xlabel('$y$ \qquad','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
axis equal;
axis([0 b_wo_wl/2 xLE_root_0-1 xLE_tip+c_tip+1])

subplot(1,3,3)
plot(y,c_equiv_vs_y)
grid minor
title('Distribuzione di corde ala equivalente','Interpreter','latex')
xlabel('$y$ \qquad','Interpreter','latex')
ylabel('$c(y)$','Interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
axis equal;
axis([0 b_wo_wl/2 0 c_root_0+1])

figure;
plot(y,gamma_a1)
grid minor
title('Carico addizionale per CL=1','Interpreter','latex')
xlabel('$y$ \qquad','Interpreter','latex')
ylabel('$\gamma$','Interpreter','latex')

figure;
plot(y,Cl_a1)
grid minor
title('Distribuzione di Cl per CL=1','Interpreter','latex')
xlabel('$y$ \qquad','Interpreter','latex')
ylabel('$Cl$','Interpreter','latex')

figure;
plot(y,Cl_a1,y,Cl_sweep)
grid minor
title('Confronto distribuzione Cl Schrenk Pope e Haney','Interpreter','latex')
xlabel('$y$ \qquad','Interpreter','latex')
ylabel('$Cl$','Interpreter','latex')

figure;
plot(y,Cl_a1,y,Cl_sweep_correct)
grid minor
title('Confronto distribuzione Cl Schrenk Pope e Haney per CL=1','Interpreter','latex')
xlabel('$y$ \qquad','Interpreter','latex')
ylabel('$Cl$','Interpreter','latex')

figure;
plot(y,Cl_b, y,Cl_a1, y, Cl_tot)
grid minor
title('Distribuzione di Cl basico e addizionale per CL=1','Interpreter','latex')
xlabel('$y$ \qquad','Interpreter','latex')
ylabel('$Cl$','Interpreter','latex')

figure;
plot(y,twist_vs_y, y,aol_vs_y, y, Cl_tot)
grid minor
title('Distribuzione di Cl basico e addizionale per CL=1','Interpreter','latex')
xlabel('$y$ \qquad','Interpreter','latex')
ylabel('$Cl$','Interpreter','latex')


c_vs_eta=[y'/(b_wo_wl/2),c_vs_y'];
c_ell_vs_eta=[y'/(b_wo_wl/2),c_ell'];
c_mean_vs_eta=[y'/(b_wo_wl/2),cmean'];
gamma_a1_vs_eta=[y'/(b_wo_wl/2),gamma_a1'];
Cl_CL_1_vs_eta=[y'/(b_wo_wl/2),Cl_a1'];
Cl_CL_04_vs_eta=[y'/(b_wo_wl/2),Cl_04'];
Cl_sweep_correct_vs_eta=[y'/(b_wo_wl/2),Cl_sweep_correct'];
gamma_sweep_vs_eta=[y'/(b_wo_wl/2),gamma_sweep_correct'];
xLE=[-flipud(y'),flipud(xLE_vs_y'-xLE_root_0);y',xLE_vs_y'-xLE_root_0];
xTE=[-flipud(y'),flipud(xLE_vs_y'-xLE_root_0+c_vs_y'); y',xLE_vs_y'-xLE_root_0+c_vs_y'];

save('c_vs_eta.dat','c_vs_eta','-ascii')
save('c_ell_vs_eta.dat','c_ell_vs_eta','-ascii')
save('c_mean_vs_eta.dat','c_mean_vs_eta','-ascii')
save('gamma_a1_vs_eta.dat','gamma_a1_vs_eta','-ascii')
save('Cl_CL_1_vs_eta.dat','Cl_CL_1_vs_eta','-ascii')
save('Cl_CL_04_vs_eta.dat','Cl_CL_04_vs_eta','-ascii')
save('Cl_sweep_correct_vs_eta.dat','Cl_sweep_correct_vs_eta','-ascii')
save('gamma_sweep_vs_eta.dat','gamma_sweep_vs_eta','-ascii')
save('xLE.dat','xLE','-ascii')
save('xTE.dat','xTE','-ascii')


function chord=c(y,c_root_0,c_kink,c_tip,yLE_kink,yLE_tip,cracked)
if cracked==1
    if y<yLE_kink
        chord=c_root_0+((c_kink-c_root_0)/(yLE_kink))*y;
    else
        chord=c_kink+((c_tip-c_kink)/(yLE_tip-yLE_kink))*(y-yLE_kink);
    end    
else 
    chord=c_root_0+((c_tip-c_root_0)/(yLE_tip))*y;
        
end
end