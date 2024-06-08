clc
clear
figure()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          30-0.4                                         %   
%                                                                         %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A2 = importdata('F:\Newexprimentstest\Trajectory\experiment_30_04.txt');
A2(:,1)=A2(:,1)-3.5;
A2(:,2)=A2(:,2)-0.4;

rou=2650;%密度2650kg/m3
d=0.000228;%平均粒径，
miu=0.0000183;%18.3μpa*s
V0=2;%2.1速度m/s
g=9.8;%重力加速度
gama=294;%图6a斜率，-300s-1
fai=34/180*pi;%图5a角度，34°
c=3.1;%图6b斜率，0.31m/s
T1=0.004;%时间，0.4ms

r1=0.01;%溅射起点距离，m
r2=0.04
r0=0.00001;%计算步长

A1=zeros(3001,3);%%%%case1(时间，x值，y值)
A1(:,1)=(r1:r0:r2);
for k=1:1:3001
    A1(k,2)=(V0-gama*A1(k,1)*sin(fai))*cos(fai)*(T1-A1(k,1)/c)-A1(k,1)+0.018;
%     A1(k,3)=-0.5*g*(T-A1(k,1)/c)^2+(V0-gama*A1(k,1)*sin(fai))*sin(fai)*(T-A1(k,1)/c);
A1(k,3)=(V0-gama*A1(k,1)*sin(fai))*sin(fai)*(T1-A1(k,1)/c);
end

A=zeros(3001,3);%%%%case(时间，x值，y值)
A(:,1)=(r1:r0:r2);
for i=1:1:3001
    A(i,2)=(rou*d^2*(V0-gama*A(i,1)*sin(fai))*cos(fai))/(18*miu)*(1-exp(-18*miu*(T1-A(i,1)/c)/rou/d^2))-A(i,1)+0.018;
    A(i,3)=((rou*d^2*(V0-gama*A(i,1)*sin(fai))*sin(fai))/18/miu+rou^2*d^4*g/324/miu^2)*(1-exp(-18*miu*(T1-A(i,1)/c)/rou/d^2))-rou*d^2*g*(T1-A(i,1)/c)/18/miu;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          30-0.8                                         %   
%                                                                         %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A3 = importdata('F:\Newexprimentstest\Trajectory\experiment_30_08.txt');
A3(:,1)=A3(:,1)-3.5;
A3(:,2)=A3(:,2)-0.4;

T2=0.008;%时间，0.4ms


r1=0.01;%溅射起点距离，m
r2=0.04
r0=0.00001;%计算步长
A4=zeros(3001,3);%%%%case(时间，x值，y值)
A4(:,1)=(r1:r0:r2);
for i=1:1:3001
    A4(i,2)=(rou*d^2*(V0-gama*A4(i,1)*sin(fai))*cos(fai))/(18*miu)*(1-exp(-18*miu*(T2-A4(i,1)/c)/rou/d^2))-A4(i,1)+0.018;
    A4(i,3)=((rou*d^2*(V0-gama*A4(i,1)*sin(fai))*sin(fai))/18/miu+rou^2*d^4*g/324/miu^2)*(1-exp(-18*miu*(T2-A4(i,1)/c)/rou/d^2))-rou*d^2*g*(T2-A(i,1)/c)/18/miu;
end

A5=zeros(3001,3);%%%%case1(时间，x值，y值)
A5(:,1)=(r1:r0:r2);
for k=1:1:3001
    A5(k,2)=(V0-gama*A5(k,1)*sin(fai))*cos(fai)*(T2-A5(k,1)/c)-A5(k,1)+0.018;
%     A5(k,3)=-0.5*g*(T-A5(k,1)/c)^2+(V0-gama*A5(k,1)*sin(fai))*sin(fai)*(T-A5(k,1)/c);
A5(k,3)=(V0-gama*A5(k,1)*sin(fai))*sin(fai)*(T2-A5(k,1)/c);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          30-1.3                                         %   
%                                                                         %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A6= importdata('F:\Newexprimentstest\Trajectory\experiment_30_13.txt');
A6(:,1)=A6(:,1)-3.5;
A6(:,2)=A6(:,2)-0.4;

T3=0.013;%时间，0.4ms

r1=0.01;%溅射起点距离，m
r2=0.04
r0=0.00001;%计算步长
A7=zeros(3001,3);%%%%case(时间，x值，y值)
A7(:,1)=(r1:r0:r2);
for i=1:1:3001
   A7(i,2)=(rou*d^2*(V0-gama*A7(i,1)*sin(fai))*cos(fai))/(18*miu)*(1-exp(-18*miu*(T3-A7(i,1)/c)/rou/d^2))-A7(i,1)+0.018;
  A7(i,3)=((rou*d^2*(V0-gama*A7(i,1)*sin(fai))*sin(fai))/18/miu+rou^2*d^4*g/324/miu^2)*(1-exp(-18*miu*(T3-A7(i,1)/c)/rou/d^2))-rou*d^2*g*(T3-A7(i,1)/c)/18/miu;
end

A8=zeros(3001,3);%%%%case1(时间，x值，y值)
A8(:,1)=(r1:r0:r2);
for k=1:1:3001
    A8(k,2)=(V0-gama*A8(k,1)*sin(fai))*cos(fai)*(T3-A8(k,1)/c)-A8(k,1)+0.018;
%     A1(k,3)=-0.5*g*(T-A1(k,1)/c)^2+(V0-gama*A1(k,1)*sin(fai))*sin(fai)*(T-A1(k,1)/c);
A8(k,3)=(V0-gama*A8(k,1)*sin(fai))*sin(fai)*(T3-A8(k,1)/c);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          30-1.9                                         %   
%                                                                         %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A9 = importdata('F:\Newexprimentstest\Trajectory\experiment_30_19.txt');
A9(:,1)=A9(:,1)-3.5;
A9(:,2)=A9(:,2)-0.4;


T4=0.019;%时间，0.4ms
r1=0.01;%溅射起点距离，m
r2=0.04
r0=0.00001;%计算步长
A10=zeros(3001,3);%%%%case(时间，x值，y值)
A10(:,1)=(r1:r0:r2);
for i=1:1:3001
    A10(i,2)=(rou*d^2*(V0-gama*A10(i,1)*sin(fai))*cos(fai))/(18*miu)*(1-exp(-18*miu*(T4-A10(i,1)/c)/rou/d^2))-A10(i,1)+0.018;
    A10(i,3)=((rou*d^2*(V0-gama*A10(i,1)*sin(fai))*sin(fai))/18/miu+rou^2*d^4*g/324/miu^2)*(1-exp(-18*miu*(T4-A10(i,1)/c)/rou/d^2))-rou*d^2*g*(T4-A10(i,1)/c)/18/miu;
end

A11=zeros(3001,3);%%%%case1(时间，x值，y值)
A11(:,1)=(r1:r0:r2);
for k=1:1:3001
    A11(k,2)=(V0-gama*A11(k,1)*sin(fai))*cos(fai)*(T4-A11(k,1)/c)-A11(k,1)+0.018;
%     A1(k,3)=-0.5*g*(T-A1(k,1)/c)^2+(V0-gama*A1(k,1)*sin(fai))*sin(fai)*(T-A1(k,1)/c);
A11(k,3)=(V0-gama*A11(k,1)*sin(fai))*sin(fai)*(T4-A11(k,1)/c);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          30-3                                         %   
%                                                                         %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A12 = importdata('F:\Newexprimentstest\Trajectory\experiment_30_30.txt');

A12 = importdata('F:\Newexprimentstest\Trajectory\experiment_30_30.txt');
A12(:,1)=A12(:,1)-3.5;
A12(:,2)=A12(:,2)-0.4;


T5=0.03;%时间，0.4ms

r1=0.01;%溅射起点距离，m
r2=0.04
r0=0.00001;%计算步长
A13=zeros(3001,3);%%%%case(时间，x值，y值)
A13(:,1)=(r1:r0:r2);
for i=1:1:3001
    A13(i,2)=(rou*d^2*(V0-gama*A13(i,1)*sin(fai))*cos(fai))/(18*miu)*(1-exp(-18*miu*(T5-A13(i,1)/c)/rou/d^2))-A13(i,1)+0.018;
    A13(i,3)=((rou*d^2*(V0-gama*A13(i,1)*sin(fai))*sin(fai))/18/miu+rou^2*d^4*g/324/miu^2)*(1-exp(-18*miu*(T5-A13(i,1)/c)/rou/d^2))-rou*d^2*g*(T5-A13(i,1)/c)/18/miu;
end

A14=zeros(3001,3);%%%%case1(时间，x值，y值)
A14(:,1)=(r1:r0:r2);
for k=1:1:3001
    A14(k,2)=(V0-gama*A14(k,1)*sin(fai))*cos(fai)*(T5-A14(k,1)/c)-A14(k,1)+0.018;
%     A1(k,3)=-0.5*g*(T-A1(k,1)/c)^2+(V0-gama*A1(k,1)*sin(fai))*sin(fai)*(T-A1(k,1)/c);
A14(k,3)=(V0-gama*A14(k,1)*sin(fai))*sin(fai)*(T5-A14(k,1)/c);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          30-3.5                                         %   
%                                                                         %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


A15= importdata('F:\Newexprimentstest\Trajectory\experiment_30_35.txt');
A15(:,1)=A15(:,1)-3.5;
A15(:,2)=A15(:,2)-0.4;


T6=0.035;%时间，0.4ms

r1=0.01;%溅射起点距离，m
r2=0.04
r0=0.00001;%计算步长
A16=zeros(3001,3);%%%%case(时间，x值，y值)
A16(:,1)=(r1:r0:r2);
for i=1:1:3001
    A16(i,2)=(rou*d^2*(V0-gama*A16(i,1)*sin(fai))*cos(fai))/(18*miu)*(1-exp(-18*miu*(T6-A16(i,1)/c)/rou/d^2))-A16(i,1)+0.018;
    A16(i,3)=((rou*d^2*(V0-gama*A16(i,1)*sin(fai))*sin(fai))/18/miu+rou^2*d^4*g/324/miu^2)*(1-exp(-18*miu*(T6-A16(i,1)/c)/rou/d^2))-rou*d^2*g*(T6-A16(i,1)/c)/18/miu;
end

A17=zeros(3001,3);%%%%case1(时间，x值，y值)
A17(:,1)=(r1:r0:r2);
for k=1:1:3001
    A17(k,2)=(V0-gama*A17(k,1)*sin(fai))*cos(fai)*(T6-A17(k,1)/c)-A17(k,1)+0.018;
%     A1(k,3)=-0.5*g*(T-A1(k,1)/c)^2+(V0-gama*A1(k,1)*sin(fai))*sin(fai)*(T-A1(k,1)/c);
A17(k,3)=(V0-gama*A17(k,1)*sin(fai))*sin(fai)*(T6-A17(k,1)/c);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% A3=-[3.39 0.89;16.41 9.66];%实验
% subplot(131)
% hold on
% plot(A(:,1).*100,A(:,2).*100,'displayname','case2','color','blue','linewidth',1);%单位，mm，m/s
% plot(A1(:,1)*100,A1(:,2)*100,'displayname','case1','color','red','linewidth',1);%单位，mm，m/s
% xlabel('R(cm)');ylabel('X(cm)');box on
% legend('show')

% subplot(132)
% hold on
% plot(A(:,1)*100,A(:,3)*100,'displayname','case2','color','blue','linewidth',1);%单位，mm，m/s
% plot(A1(:,1)*100,A1(:,3)*100,'displayname','case1','color','red','linewidth',1);%单位，mm，m/s
% xlabel('R(cm)');ylabel('Y(cm)');box on
% legend('show')
% 
% subplot(133)
% hold on

xstring='x_\psi(cm)'

ystring='y_\psi(cm)'
%  ystring2='Dn(mm)'
%  ystring3='Db(mm)'
% ystring4='H(mm)'
xstring2='x(cm)'

 ystring2='y(cm)'


figure(1);
txtsize = 24;

hold on
%     fig_plot1 = scatter(Modulus1_x,Modulus1_y);
%     fig_plot2 = scatter(Modulus2_x,Modulus2_y);
%     fig_plot3 = scatter(Modulus3_x,Modulus3_y);Anglevariation1_x
%     fig_plot4 = scatter(Modulus4_x,Modulus4_y);
%     fig_plot5 = scatter(Modulus5_x,Modulus5_y);
%     fig_plot6 = scatter(Modulus6_x,Modulus6_y);
%     fig_plot7 = scatter(Modulus7_x,Modulus7_y);
%     fig_plot8 = scatter(Modulus8_x,Modulus8_y);
    
%     fig_plot1 = scatter(log10(Modulus1_x),log10(Modulus1_y));
%     fig_plot2 = scatter(log10(Modulus2_x),log10(Modulus2_y));
%     fig_plot3 = scatter(log10(Modulus3_x),log10(Modulus3_y));
%     fig_plot4 = scatter(log10(Modulus4_x),log10(Modulus4_y));
%     fig_plot5 = scatter(log10(Modulus5_x),log10(Modulus5_y));
% %     fig_plot6 = scatter(log10(Modulus6_x),log10(Modulus6_y));
% %     fig_plot7 = scatter(log10(Modulus7_x),log10(Modulus7_y));
%     fig_plot8 = scatter(log10(Modulus8_x),log10(Modulus8_y));
% 

%  fig_plot1 =scatter(Modulus8_x,Modulus8_y)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot1 = figure('OuterPosition',[400,50,1500,600]);

     subplot(1,2,1)
     fig_plot1 =scatter(A2(:,1),A2(:,2))
     hold on

     fig_plot2 =scatter(A3(:,1),A3(:,2));
     hold on
    
    fig_plot3 =scatter(A6(:,1),A6(:,2));
    hold on
    fig_plot4 =scatter(A9(:,1),A9(:,2));
    hold on
    fig_plot5 = scatter(A12(:,1),A12(:,2));
    hold on
    fig_plot6 = scatter(A15(:,1),A15(:,2));
    
      hold on

    fig_plot7 =plot(A1(:,2)*100,A1(:,3)*100,'--','LineWidth',1,'Color',[0 0 0])
    hold on

    fig_plot8= plot(A5(:,2)*100,A5(:,3)*100,'--','LineWidth',1,'Color',[0 0 0]);
    hold on
    
    fig_plot9 = plot(A8(:,2)*100,A8(:,3)*100,'--','LineWidth',1,'Color',[0 0 0]);
    hold on
    fig_plot10 = plot(A11(:,2)*100,A11(:,3)*100,'--','LineWidth',1,'Color',[0 0 0]);
    hold on
    fig_plot11 = plot(A14(:,2)*100,A14(:,3)*100,'--','LineWidth',1,'Color',[0 0 0]);
    hold on
    
    fig_plot12 = plot(A17(:,2)*100,A17(:,3)*100,'--','LineWidth',1,'Color',[0 0 0]);
    hold on
    
    fig_plot13 = plot(A(:,2)*100,A(:,3)*100,'-','LineWidth',1,'Color',[0 0 0]);
    hold on
        
    fig_plot14 = plot(A4(:,2)*100,A4(:,3)*100,'-','LineWidth',1,'Color',[0 0 0]);
    hold on
        
    fig_plot15 = plot(A7(:,2)*100,A7(:,3)*100,'-','LineWidth',1,'Color',[0 0 0]);
    hold on
        
    fig_plot16 = plot(A10(:,2)*100,A10(:,3)*100,'-','LineWidth',1,'Color',[0 0 0]);
    hold on
        
    fig_plot17 = plot(A13(:,2)*100,A13(:,3)*100,'-','LineWidth',1,'Color',[0 0 0]);
    hold on   
    fig_plot18 = plot(A16(:,2)*100,A16(:,3)*100,'-','LineWidth',1,'Color',[0 0 0]);
    hold on
    
    
    
    
     set(gca,'linewidth',2, 'FontSize',15, 'Box', 'on');
%      axis([0,18,-2,9.8])%up
%      set(gca,'XTick',[]);
           axis([-0.2,30,0.4,30])%neck
                    
%            axis([0,4.5,1,16])%bottom
%            axis([0,4.5,0,9])%Height
    legend('4ms','8ms','13ms','19ms','30ms','35ms','Location','SouthEast','FontSize',15)
%     legend('(a?','Location','SouthEast','FontSize',15)
    set( legend,'Box', 'off')
    xlabel(xstring, 'FontName','Times New Roman', 'FontSize', txtsize)
    ylabel(ystring, 'FontName','Times New Roman', 'FontSize', txtsize)
    set(fig_plot1,'LineWidth',3,'Marker','+')%'MarkerSize',12)
    set(fig_plot2,'LineWidth',3,'Marker','o')%'MarkerSize',12)
    set(fig_plot3,'LineWidth',3,'Marker','s')%'MarkerSize',12)
    set(fig_plot4,'LineWidth',3,'Marker','v')%'MarkerSize',12)
    set(fig_plot5,'LineWidth',3,'Marker','p')%'MarkerSize',12)
    set(fig_plot6,'LineWidth',3,'Marker','h')%'MarkerSize',12)
    legend('(a)','FontSize',15)

    hold on
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     
%     

