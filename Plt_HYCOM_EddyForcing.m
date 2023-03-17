format compact, clear

Tmod=[0:1/24:7-1/24]; L=length(Tmod);
load vertical_grid_128
z=vertical_grid_128(2:end-1,2);

load LSF_TSEquation_C.dat
T_Eddy=LSF_TSEquation_C(:,3); S_Eddy=LSF_TSEquation_C(:,6);
T2D=reshape(T_Eddy,128,L); S2D=reshape(S_Eddy,128,L);

load LSF_MomentumEquation_C.dat
U_Eddy=LSF_MomentumEquation_C(:,5); V_Eddy=LSF_MomentumEquation_C(:,9);
UV=sqrt(U_Eddy.^2+V_Eddy.^2);
UV2D=reshape(UV,128,L);

colormap('jet');
subplot(3,1,1)
contourf(Tmod,z,T2D*1.0e6,40,'linecolor','none');
h1=colorbar;
title(h1,'^oC s^-^1 x 10^-^6');
hold on
plot([62/24 62/24],[-200 0],'w-','linewidth',1);
plot([65/24 65/24],[-200 0],'w-','linewidth',1);
title('Temperature Eddy Forcing at Point C','fontsize',20);
set(gca,'xlim',[1/2 5],'xtick',[1/2:1/2:5],'xticklabel',[12:12:120]);
set(gca,'ylim',[-80 0],'ytick',[-80:10:0],'yticklabel',[-80:10:0]);
set(gca,'fontsize',16);
ylabel('Depth (m)');

subplot(3,1,2)
contourf(Tmod,z,S2D*1.0e6,40,'linecolor','none');
h2=colorbar;
title(h2,'psu s^-^1 x 10^-^6');
hold on
plot([62/24 62/24],[-200 0],'w-','linewidth',1);
plot([65/24 65/24],[-200 0],'w-','linewidth',1);
title('Salinity Eddy Forcing at Point C','fontsize',20);
set(gca,'xlim',[1/2 5],'xtick',[1/2:1/2:5],'xticklabel',[12:12:120]);
set(gca,'ylim',[-80 0],'ytick',[-80:10:0],'yticklabel',[-80:10:0]);
set(gca,'fontsize',16);
ylabel('Depth (m)');

subplot(3,1,3)
contourf(Tmod,z,UV2D*1.0e8,40,'linecolor','none');
h3=colorbar;
title(h3,'m s^-^2 x 10^-^8');
hold on
plot([62/24 62/24],[-200 0],'w-','linewidth',1);
plot([65/24 65/24],[-200 0],'w-','linewidth',1);
title('Momentum Eddy Forcing at Point C','fontsize',20);
set(gca,'xlim',[1/2 5],'xtick',[1/2:1/2:5],'xticklabel',[12:12:120]);
set(gca,'ylim',[-80 0],'ytick',[-80:10:0],'yticklabel',[-80:10:0]);
set(gca,'fontsize',16);
ylabel('Depth (m)');
xlabel('Time (hours)');


