clc; close all; clear all;
format long

%Dati velivolo
b = 58.30 %m
S = 370.45 %m
AR = b^2/S 

c_root = 11.98;    xLE_root=18.626;      yLE_root=0.0;
c_kink = 7.200;    xLE_kink=24.35;        yLE_kink=9.1; 
c_tip = 2.57;      xLE_tip=36.973;        yLE_tip=29.15; 

MLW = 267000 %KG
OEW = 129500 %KG 
MTOW = 275000 %kg

n_e = 4; %numero di motori sotto l'ala
V_min = 247 / 3.6; %m/s
V_cruise = 630/3.6 %m/s

c_tip = 2.57 %m
Sweep_LE = atan((xLE_tip-xLE_root)/(b/2)) *(180)/pi; %deg
Sweep_LE_rad = atan((xLE_tip-xLE_root)/(b/2));


%Dati condizioni operative
g= 9.81 %m/s^2
rho_sl = 1.225 %kg/m^3
T0= 15 %C
k_exp = 4.256
Tz= 0.0065 %°/m
h_vec = 0:100:11000;

%Derived data
y=linspace(0,b/2,100);
eta=linspace(0,1,100);
c_vs_y=nan(1,length(y));
c_equiv_vs_y=nan(1,length(y));
for i=1:length(y)
c_vs_y(i)=c(y(i),c_root,c_kink,c_tip,yLE_kink,yLE_tip,1);
end
xLE_vs_y=xLE_root+((xLE_tip-xLE_root)/(yLE_tip-yLE_root))*(y-yLE_root);
xTE_vs_y=xLE_vs_y+c_vs_y;
y_wing_planform = [ yLE_root; yLE_kink; yLE_tip; yLE_tip; yLE_kink; yLE_root];
x_wing_planform = [0;xLE_kink - xLE_root; xLE_tip - xLE_root; xLE_tip - xLE_root + c_tip; xLE_kink - xLE_root + c_kink; c_root];

%-----------------------------------------------------------------------
%ALA EQUIVALENTE
c_root_equiv=S*2/b-c_tip;
taperRatio_equiv=c_tip/c_root_equiv;
xTE_equiv_vs_y=xLE_vs_y+c_equiv_vs_y;
for i=1:length(y)
c_equiv_vs_y(i)=c(y(i),c_root_equiv,c_kink,c_tip,yLE_kink,yLE_tip,0);
end
y_wing_planform_eq = [yLE_root; yLE_tip; yLE_tip; yLE_root];
x_wing_planform_eq = [0; xLE_tip - xLE_root; xLE_tip - xLE_root + c_tip; c_root_equiv];

%CL MAX
CL_max_lan_sl = (MLW  * g)/(0.5* rho_sl * V_min^2 *S);

%V min
%V min al variare della quota a MLW

Vmin_vs_h =nan(1,length(h_vec));
for i=1:length(h_vec)
    rho = rhoCalc(h_vec(i),k_exp, T0, Tz, rho_sl);
    Vmin_vs_h(i)= sqrt((MLW  * g)/(0.5* rho * CL_max_lan_sl *S)); %NB. dovrebbe essere iterativo
end


% CL
CL_vs_h_MLW =nan(1,length(h_vec));
for i=1:length(h_vec)
    rho = rhoCalc(h_vec(i),k_exp, T0, Tz, rho_sl);
    CL_vs_h_MLW(i)=  (MLW  * g)/(0.5* rho * V_cruise^2 *S);
end

CL_vs_h_MTOW =nan(1,length(h_vec));
for i=1:length(h_vec)
    rho = rhoCalc(h_vec(i),k_exp, T0, Tz, rho_sl);
    CL_vs_h_MTOW(i)=  (MTOW  * g)/(0.5* rho * V_cruise^2 *S) ;
end

CL_vs_h_OEW =nan(1,length(h_vec));
for i=1:length(h_vec)
    rho = rhoCalc(h_vec(i),k_exp, T0, Tz, rho_sl);
    CL_vs_h_OEW(i)=  (OEW  * g)/(0.5* rho * V_cruise^2 *S) ;
end

%PLOT
plot(y_wing_planform,x_wing_planform, y_wing_planform_eq ,x_wing_planform_eq,'r')
axis ij
grid minor
title('Pianta semiala','Interpreter','latex')
xlabel('$y$ \qquad','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
axis equal;

figure
plot(h_vec,Vmin_vs_h)
grid minor
title('Velocità minima al variare della quota','Interpreter','latex')
xlabel('$h$ \qquad','Interpreter','latex')
ylabel('$CL$','Interpreter','latex')
set(get(gca,'ylabel'))

figure
plot(h_vec,CL_vs_h_MLW,h_vec,CL_vs_h_MTOW,h_vec,CL_vs_h_OEW)
grid minor
title('CL al variare della quota e del peso','Interpreter','latex')
xlabel('$h$ \qquad','Interpreter','latex')
ylabel('$C_L$','Interpreter','latex')
set(get(gca,'ylabel'))




%Outputs-----
wing_planform=[y_wing_planform',x_wing_planform'];
wing_planform_eq=[y_wing_planform_eq',x_wing_planform_eq'];
vmin_vs_h = [h_vec',Vmin_vs_h'];
cl_vs_h_OEW =[h_vec', CL_vs_h_OEW'];
cl_vs_h_MLW =[h_vec', CL_vs_h_MLW'];
cl_vs_h_MTOW =[h_vec', CL_vs_h_MTOW'];


save('1_wing_planform.dat','wing_planform','-ascii')
save('1_wing_planform_eq.dat','wing_planform_eq','-ascii')
save('2_Vmin_vs_h_CLmaxL_fixed.dat','vmin_vs_h','-ascii')
save('3_Cl_vs_h_fixedW_OEW.dat','cl_vs_h_OEW','-ascii')
save('3_Cl_vs_h_fixedW_MLW.dat','cl_vs_h_MLW','-ascii')
save('3_Cl_vs_h_fixedW_MTOW.dat','cl_vs_h_MTOW','-ascii')


%Functions----
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

function rho = rhoCalc(h,k_exp, T0, Tz, rho0)
T_z = (T0 - Tz*h)+ 273.15;
rho = rho0*(T_z/(T0 + 273.15))^(k_exp)
end