clc
clear
figure()


I7 = importdata('D:\2021-Ejecta experiment\Propogationvelocity.txt');




Propogationvelocity1_x = I7(:,1)+0.14; Propogationvelocity1_y = I7(:,2);
Propogationvelocity2_x = I7(:,3)+0.2; Propogationvelocity2_y = I7(:,4);
Propogationvelocity3_x = I7(:,5)+0.8; Propogationvelocity3_y = I7(:,6);
Propogationvelocity4_x = I7(:,7)+0.7; Propogationvelocity4_y = I7(:,8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          30                                             %   
%                                                                         %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rou=2650;%密度2650kg/m3
d=0.000228;%平均粒径，
miu=0.0000183;%18.3μpa*s
V0=2.1;%2.1速度m/s
g=9.8;%重力加速度
h=0.0001;
gama=30.4;%图6a斜率，-300s-1
fai=29/180*pi;%图5a角度，34°
c=5.5;%图6b斜率，0.31m/s


r1=0.0001;
r2=0.1428;
r0=0.001;%计算步长
N=fix((r2*100000-r1*100000)/(r0*100000))+1
T=zeros(N,2);
T(:,1)=(r1:r0:r2);
for k=1:1:N
    T(k,2)=T(k,1)/c+h/(sin(fai)*V0*(1-T(k,1)*gama*sin(fai)/V0)); 
end


x1=29.76;x2=60;y1=141.1;
M=fix((x2*10000-x1*10000)/(r0*10000))+1;
sita=11.3/180*pi;
hh=10;
kk=0.001;
gamaa=4;
[t,r]=ode45(@(t,r) gamaa*hh*(log(r/hh)-sita*r/hh-kk),(x1:r0:x2),y1); 
% plot(t,r,'k-') 
% hold on

TR1=zeros(M+N,2);
TR1(:,1)=[T(:,1);r];
TR1(:,2)=[T(:,2);t];

 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                             45                                          %   
% %                                                                         %   
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rou1=2650;%密度2650kg/m3
d1=0.000228;%平均粒径，
miu1=0.0000183;%18.3μpa*s
V01=3.0;%2.1速度m/s
g1=9.8;%重力加速度
h1=0.001;
gama1=38;%图6a斜率，-300s-1
fai1=45/180*pi;%图5a角度，34°
c1=5;%图6b斜率，0.31m/s


r3=0.0001;
r4=0.109;
r0=0.0001;%计算步长
N1=fix((r4*1000-r3*1000)/(r0*1000))+1
T1=zeros(N1,2);
T1(:,1)=(r3:r0:r4);
for k=1:1:N1
    T1(k,2)=T1(k,1)/c1+h1/(sin(fai1)*V01*(1-T1(k,1)*gama1*sin(fai1)/V01));
end


x12=31;x22=60;y12=108.24;
M1=fix((x22*10000-x12*10000)/(r0*10000))+1;
sita=13/180*pi;
hh=9;
k=0.001;
[t,r]=ode45(@(t,r) gamaa*hh*(log(r/hh)-sita*r/hh-k),(x12:r0:x22),y12); 
% plot(t,r,'k-') 
%     hold on

TR2=zeros(M1+N1,2);
TR2(:,1)=[T1(:,1);r];
TR2(:,2)=[T1(:,2);t];



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                          60                                            %   
% %                                                                         %   
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
rou2=2650;%密度2650kg/m3
d2=0.000228;%平均粒径，
miu2=0.0000183;%18.3μpa*s
V02=1.8;%2.1速度m/s
g2=9.8;%重力加速度
h2=0.0005;
gama2=22;%图6a斜率，-300s-1
fai2=52/180*pi;%图5a角度，34°
c2=4;%图6b斜率，0.31m/s


r5=0.0001;
r6=0.103;
r0=0.001;%计算步长
N2=fix((r6*1000-r5*1000)/(r0*1000))+1
T2=zeros(N2,2);
T2(:,1)=(r5:r0:r6);
for k=1:1:N2
    T2(k,2)=T2(k,1)/c2+h2/(sin(fai2)*V02*(1-T2(k,1)*gama2*sin(fai2)/V02));
end

x13=30.9;x23=60;y13=98;
M2=fix((x23*10000-x13*10000)/(r0*10000))+1;
sita=15/180*pi;
hh=10;
k=0.001;
[t,r]=ode45(@(t,r) gamaa*hh*(log(r/hh)-sita*r/hh-k),(x13:r0:x23),y13); 
% plot(t,r,'k-') 
%     hold on

TR3=zeros(M2+N2,2);
TR3(:,1)=[T2(:,1);r];
TR3(:,2)=[T2(:,2);t];

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                          90                                            %   
% %                                                                         %   
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
rou3=2650;%密度2650kg/m3
d3=0.000228;%平均粒径，
miu3=0.0000183;%18.3μpa*s
V03=1.1;%2.1速度m/s
g3=9.8;%重力加速度
h3=0.0001;
gama3=19.7;%图6a斜率，-300s-1
fai3=52.5/180*pi;%图5a角度，34°
c3=3;%图6b斜率，0.31m/s


r7=0.0001;
r8=0.07108;
r0=0.001;%计算步长
N3=fix((r8*100000-r7*100000)/(r0*100000))+1
T3=zeros(N3,2);
T3(:,1)=(r7:r0:r8);
for k=1:1:N3
    T3(k,2)=T3(k,1)/c3+h3/(sin(fai3)*V03*(1-T3(k,1)*gama3*sin(fai3)/V03));
end

x14=26;x24=60;y14=68;
M3=fix((x24*10000-x14*10000)/(r0*10000))+1;
sita=15.5/180*pi;
hh=7;
k=0.001;
[t,r]=ode45(@(t,r) gamaa*hh*(log(r/hh)-sita*r/hh-k),(x14:r0:x24),y14); 
% plot(t,r,'k-') 
%     hold on

TR4=zeros(M3+N3,2);
TR4(:,1)=[T3(:,1);r];
TR4(:,2)=[T3(:,2);t];
% 
% 
% 
% plot1 = figure('OuterPosition',[400,50,1500,600]);

% fig_plot1 =plot(T(:,2)*1000,T(:,1)*1000,'--','LineWidth',1,'Color',[0 0 0]);
% hold on
% % 
% % fig_plot2 =plot(T1(:,2)*1000,T1(:,1)*1000,'--','LineWidth',1,'Color',[0 0 0]);
% % hold on
% % 
% % fig_plot3 =plot(T2(:,2)*1000,T2(:,1)*1000,'--','LineWidth',1,'Color',[0 0 0]);
% % hold on
% % 
% % fig_plot4 =plot(T3(:,2)*1000,T3(:,1)*1000,'--','LineWidth',1,'Color',[0 0 0]);
% % hold on
% % set(gca,'linewidth',2, 'FontSize',15, 'Box', 'on');
% % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xstring='t(ms)'
ystring='r(cm)'

figure(1);
txtsize = 24;

hold on

plot1 = figure('OuterPosition',[400,50,1500,600]);

     subplot(1,2,1)
plot(T(:,2)*1000,T(:,1)*100,'--s','LineWidth',1,'MarkerSize',3,'Color','k');
hold on

fig_plot2 =plot(T1(:,2)*1000,T1(:,1)*100,'b--o','LineWidth',1,'MarkerSize',3,'Color','c');
hold on

fig_plot3 =plot(T2(:,2)*1000,T2(:,1)*100,'b--*','LineWidth',1,'Color','r');
hold on

fig_plot4 =plot(T3(:,2)*1000,T3(:,1)*100,'--gs','LineWidth',1,'MarkerSize',3,'Color','g');
hold on
set(gca,'linewidth',2, 'FontSize',15, 'Box', 'on');


 set( legend,'Box', 'off')
     xlabel(xstring, 'FontName','Times New Roman', 'FontSize', txtsize)
     ylabel(ystring, 'FontName','Times New Roman', 'FontSize', txtsize)
   
hold on
     
    fig_plot5 =scatter(Propogationvelocity1_x*10,Propogationvelocity1_y)
     hold on

     fig_plot6 = scatter(Propogationvelocity2_x*10,Propogationvelocity2_y);
     hold on
    
    fig_plot7 = scatter(Propogationvelocity3_x*10,Propogationvelocity3_y);
    hold on
    fig_plot8 = scatter(Propogationvelocity4_x*10,Propogationvelocity4_y);
    hold on

    set(fig_plot5,'LineWidth',3,'Marker', '+' )%'MarkerSize',12)
    set(fig_plot6,'LineWidth',3,'Marker','o')%'MarkerSize',12)
    set(fig_plot7,'LineWidth',3,'Marker','s')%'MarkerSize',12)
    set(fig_plot8,'LineWidth',3,'Marker','v')%'MarkerSize',12)

    
    
   hold on
set(gca,'linewidth',2, 'FontSize',15, 'Box', 'on');

legend('Theory for  30^\circ',' 45^\circ', '60^\circ','90^\circ','Experiment for 30^\circ','45^\circ', '60^\circ','90^\circ','FontSize',15)

 set( legend,'Box', 'off')
     xlabel(xstring, 'FontName','Times New Roman', 'FontSize', txtsize)
     ylabel(ystring, 'FontName','Times New Roman', 'FontSize', txtsize)
     
     axis([0,50,0,15])
     
% plot2 = figure('OuterPosition',[400,50,1500,600]);
subplot(1,2,2)

txtsize = 24;
plot(TR1(:,2),TR1(:,1)/10,'--s','LineWidth',1,'MarkerSize',3,'Color','k');
hold on

fig_plot2 =plot(TR2(:,2),TR2(:,1)/10,'b--o','LineWidth',1,'MarkerSize',3,'Color','c');
hold on

fig_plot3 =plot(TR3(:,2),TR3(:,1)/10,'b--*','LineWidth',1,'Color','r');
hold on

fig_plot4 =plot(TR4(:,2),TR4(:,1)/10,'--gs','LineWidth',1,'MarkerSize',3,'Color','g');
hold on
set(gca,'linewidth',2, 'FontSize',15, 'Box', 'on');

% legend('30^\circ','45^\circ', '60^\circ','90^\circ','FontSize',15)

%  set( legend,'Box', 'off')
     xlabel(xstring, 'FontName','Times New Roman', 'FontSize', txtsize)
     ylabel(ystring, 'FontName','Times New Roman', 'FontSize', txtsize)


     
% plot2 = figure('OuterPosition',[400,50,1500,600]);
% 
% fig_plot5 =plot(TR1(:,2),TR1(:,1),'--','LineWidth',1,'Color','k');
% hold on
% % 
% fig_plot6 =plot(TR2(:,2),TR2(:,1),'--','LineWidth',1,'Color','r');
% hold on
% % 
% fig_plot7 =plot(TR3(:,2)*1000,TR3(:,1),'--','LineWidth',1,'Color','g');
% hold on
% 
% fig_plot8 =plot(TR4(:,2)*1000,TR4(:,1),'--','LineWidth',1,'Color','y');
% hold on
% set(gca,'linewidth',2, 'FontSize',15, 'Box', 'on');

axis([10,50,0,15])

