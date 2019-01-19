clc; close all; clear all;

S_w=361.6;
b_tot=60.3;
b_wo_wl=58.3;
root_kinkAF=load('sc(2)-0614chiuso.dat');
tipkAF=load('sc(2)-0610chiuso.dat');

c_root_0=11.982; xLE_root_0=18.626; yLE_root_0=0.0;   zLE_root_0=-1.464;
c_root_f=10.500; xLE_root_f=20.400; yLE_root_f=2.820; zLE_root_f=-1.15;
c_kink=7.200;      xLE_kink=24.35;    yLE_kink=9.1;     zLE_kink=-0.45;
c_tip=2.57;        xLE_tip=36.973;    yLE_tip=29.15;     zLE_tip=0.86;
c_wl_tip=0.93;   xLE_wl_tip=39.427; yLE_wl_tip=30.15;  zLE_wl_tip=2.161;

c0=4*S_w/(pi*b_wo_wl);

y=linspace(0,b_wo_wl/2,100000);
cell=c0*(2/b_wo_wl)*sqrt(0.25*b_wo_wl^2-y.^2);

A1=2*trapz(y,cell)

c_vs_y=nan(1,length(y));
xTE_vs_y=nan(1,length(y));
xLE_vs_y=nan(1,length(y));
airfoil_1=nan(size(root_kinkAF,1),2,length(y));
airfoil_2=nan(size(tipkAF,1),2,length(y));
 
for i=1:length(y)
c_vs_y(i)=c(y(i),c_root_0,c_kink,c_tip,yLE_kink,yLE_tip,1);
xLE_vs_y(i)=LE(y(i),yLE_tip,xLE_root_f,xLE_tip,yLE_root_f);
xTE_vs_y(i)=LE(y(i),yLE_tip,xLE_root_f,xLE_tip,yLE_root_f)+c(y(i),c_root_0,c_kink,c_tip,yLE_kink,yLE_tip,1);
 if y(i)<yLE_kink
    airfoil_1(:,:,i)=root_kinkAF*c_vs_y(i);
 else
    airfoil_2(:,:,i)=tipkAF*c_vs_y(i); 
 end
end

A2=2*trapz(y,c_vs_y)

S_recalc=2*integral(@(y)c(y,c_root_0,c_kink,c_tip,yLE_kink,yLE_tip,1),0,b_wo_wl/2)

mac=2*integral(@(y)c2(y,c_root_0,c_kink,c_tip,yLE_kink,yLE_tip,1),0,b_wo_wl/2)/S_w

xLE_mac=2*integral(@(y)cxLE(y,c_root_0,c_kink,c_tip,yLE_kink,yLE_tip,xLE_root_f,xLE_tip,yLE_root_f,1),0,b_wo_wl/2)/S_w

yLE_mac=2*integral(@(y)cyLE(y,c_root_0,c_kink,c_tip,yLE_kink,yLE_tip,1),0,b_wo_wl/2)/S_w


figure;
plot(y,c_vs_y,y,cell,y,0.5*(c_vs_y+cell),[yLE_mac yLE_mac],[0 mac],'r')
grid minor
title('Distribuzione di corde','Interpreter','latex')
xlabel('$y$ \qquad','Interpreter','latex')
ylabel('$c(y)$','Interpreter','latex')
legend('c','c ell', '\gamma_1')
set(get(gca,'ylabel'),'rotation',0)
axis equal;
axis([0 b_wo_wl/2 0 c_root_0+1])

figure;
plot(y,xTE_vs_y,y,xLE_vs_y,[yLE_mac yLE_mac],[xLE_mac xLE_mac+mac],'r')
axis ij
grid minor
title('Pianta semiala','Interpreter','latex')
xlabel('$y$ \qquad','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
axis equal;
axis([0 b_wo_wl/2 xLE_root_0-1 xLE_tip+c_tip+1])

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

function chord2=c2(y,c_root_0,c_kink,c_tip,yLE_kink,yLE_tip,cracked)
if cracked==1
    if y<yLE_kink
        chord=c_root_0+((c_kink-c_root_0)/(yLE_kink))*y;
    else
        chord=c_kink+((c_tip-c_kink)/(yLE_tip-yLE_kink))*(y-yLE_kink);
    end    
else 
    chord=c_root_0+((c_tip-c_root_0)/(yLE_tip))*y;
        
end
chord2=chord.^2;
end     

function chord_xLE=cxLE(y,c_root_0,c_kink,c_tip,yLE_kink,yLE_tip,xLE_root_f,xLE_tip,yLE_root_f,cracked)
if cracked==1
    if y<yLE_kink
        chord=c_root_0+((c_kink-c_root_0)/(yLE_kink))*y;
    else
        chord=c_kink+((c_tip-c_kink)/(yLE_tip-yLE_kink))*(y-yLE_kink);    
    end    
else 
    chord=c_root_0+((c_tip-c_root_0)/(yLE_tip))*y;
        
end
xLE_vs_y=xLE_root_f+((xLE_tip-xLE_root_f)/(yLE_tip-yLE_root_f))*(y-yLE_root_f);
chord_xLE=chord.*xLE_vs_y;
end

function xLE=LE(y,yLE_tip,xLE_root_f,xLE_tip,yLE_root_f)

   xLE=xLE_root_f+((xLE_tip-xLE_root_f)/(yLE_tip-yLE_root_f))*(y-yLE_root_f);

end

function chord_yLE=cyLE(y,c_root_0,c_kink,c_tip,yLE_kink,yLE_tip,cracked)
if cracked==1
    if y<yLE_kink
        chord=c_root_0+((c_kink-c_root_0)/(yLE_kink))*y;
    else
        chord=c_kink+((c_tip-c_kink)/(yLE_tip-yLE_kink))*(y-yLE_kink);
    end    
else 
    chord=c_root_0+((c_tip-c_root_0)/(yLE_tip))*y;
        
end
yLE_vs_y=y;
chord_yLE=chord.*yLE_vs_y;
end