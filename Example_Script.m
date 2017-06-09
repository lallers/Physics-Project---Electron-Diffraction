%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Numerical Computation, 2016                        /
%     University of New Mexico                           /
%NOTE: None of my scripts are built to be robust, they   /
%      are merely an implementation of a given set of    /
%      data or instructions!                             /
%/////////////////////////////////////////////////////////
%Useful documents https://www.3bscientific.com/product-manual/1013885_EN.pdf

clc; clear all;close all
load('big_inner.mat');load('big_outer.mat');load('small_inner.mat');load('small_outer.mat')

big_outer = unnamed1; big_inner = unnamed;
clear unnamed unnamed1
%/////////// Constants ///////////////////////////////////
const.q = 1.6*10^(-19); %Coulombs, Charge of an Electron
const.M = 9.11*10^(-31); %kg, Electron Mass
const.h = 6.62307004*10^(-34); %Planks constant
const.hbar = 1.0545718*10^(-34);%Planks constant
const.L = 13.2; %Length of Glass
const.c = 3*10^8; %Speed of light
const.d_small_actual = 0.213; %nm
const.d_big_actual = 0.123; %nm
d =@(D,v) 4.*pi.*const.L.*const.hbar.*const.c./(D.*sqrt(2.*const.q.*v.*const.M.*const.c^2));
%///////////////////////////////////////////////////////
V = 1000*[5,4.8,4.6,4.4,4.2,4,3.5]';

%//////Don't know if needed
%theta_small_p = rad2deg(atan(([Run1(:,2),Run2(:,2),Run3(:,2),Run4(:,2)])./(2*Lp))); %Positive
%theta_small_n = rad2deg(atan(([Run1(:,2),Run2(:,2),Run3(:,2),Run4(:,2)])./(2*Ln))); %Negative
%theta_big_p = rad2deg(atan(([Run1(:,1),Run2(:,1),Run3(:,1),Run4(:,1)])./(2*Lp))); %Positive
%theta_big_n = rad2deg(atan(([Run1(:,1),Run2(:,1),Run3(:,1),Run4(:,1)])./(2*Ln))); %Negative
%///////////

data_small_D = (small_outer + small_inner)./2; %Finds center of the 'small' ring.
data_small_d = d(data_small_D,V).*10^10;
data_big_D = (big_outer + big_inner)./2; %Finds center of the 'big' ring.
data_big_d   = d(data_big_D,V).* 10^10; %Convert to small d

%////////// Do analysis to Small ring //////////////////////
mean_small_D = mean(data_small_D,2);
std_small_D = std(data_small_D')';
fit_small_D = polyfit(V,mean_small_D,1);
fit_val_small_D = polyval(fit_small_D,V);

mean_small_d = mean(data_small_d,2);
std_small_d = std(data_small_d')';
fit_small_d = polyfit(V,mean_small_d,1);
fit_val_small_d = polyval(fit_small_d,V);

final_d_small_mu = mean(fit_val_small_d);
final_d_small_sigma = std(fit_val_small_d); 
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%////////// Do analysis to Big ring //////////////////////
mean_big_D = mean(data_big_D,2);
std_big_D = std(data_big_D')';
fit_big_D = polyfit(V,mean_big_D,1);
fit_val_big_D = polyval(fit_big_D,V);

mean_big_d = mean(data_big_d,2);
std_big_d  = std(data_big_d')';
fit_big_d = polyfit(V,mean_big_d,1);
fit_val_big_d = polyval(fit_big_d,V);

final_d_big_mu = mean(fit_val_big_d);
final_d_big_sigma = std(fit_val_big_d);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%changing voltage for plotting
Vo = 1./sqrt(V);

%Plots of D v V(-1/2)
figure('Name','Diiffraction Measured Data','units','normalized','position',[.1 .1 .7 .7]);
subplot(2,2,1)
title('Small Ring (D Measured)')
hold on
plot(Vo,mean_small_D,'bo');
plot(Vo,fit_val_small_D,'r--');
errorbar(Vo,mean_small_D,std_small_D,'bo')
str = sprintf('y = %.03fx + %.03f',fit_small_D);
legend({'Mean',str},'Box','off','Location','southeast');xlabel('Voltage^{-1/2} (V)');ylabel('Diameter (mm)');
hold off

subplot(2,2,2)
title('Small Ring (d Calculated)');
hold on
plot(Vo,mean_small_d, 'bo')
plot(Vo,fit_val_small_d,'r--')
%line([min(Vo) max(Vo)],[d_small_actual, d_small_actual],'Color','m')%Actual Value
errorbar(Vo,mean_small_d,std_small_d,'o')
str = sprintf('y = %fx + %.03f',fit_small_d);
legend({'Mean',str},'Box','off','Location','southeast');xlabel('Voltage^{-1/2} (V)');ylabel('Diameter (nm)');
hold off

subplot(2,2,3)
hold on
title('Big Ring (D Measured)')
plot(Vo,mean_big_D,'bo');
plot(Vo,fit_val_big_D,'r--');
errorbar(Vo,mean_big_D,std_big_D,'bo')
str = sprintf('y = %.03fx + %.03f',fit_big_D);
legend({'Mean',str},'Box','off','Location','southeast');xlabel('Voltage^{-1/2} (V)');ylabel('Diameter (mm)');
hold off

subplot(2,2,4)
title('Big Ring (d Calculated)');
hold on
plot(Vo,mean_big_d, 'bo')
plot(Vo,fit_val_big_d,'r--')
%line([min(Vo) max(Vo)],[d_big_actual, d_big_actual],'Color','m') %Actual Value
errorbar(Vo,mean_big_d,std_big_d,'o')
str = sprintf('y = %fx + %.03f',fit_big_d);
legend({'Mean',str},'Box','off','Location','southeast');xlabel('Voltage^{-1/2} (V)');ylabel('Diameter (nm)');
hold off

%$$$$$$$$$$$$$$$$$$$$$$$$ EXTRAS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%DeBroglie Relation: lambda = h/p
velocity = sqrt(2.*const.q.*V./const.M);
wavelength = (const.h./(const.M*velocity))*10^10; %In nanometers

%/////////////// Error anaylsis stuff
% alpha = 0.05; 
% mu = mean(mean_small); %mean
% sigma = std(mean_small);   %std

% [mu, sigma, xci, yci] = normfit(mean(data_small_d));
% x = linspace(mu - 5*sigma, mu + 5*sigma,1000);
% y = normpdf(x, mu, sigma);
% 
% 
% figure('Name','Small ring calulated d','units','normalized','position',[.1 .1 .5 .7]);
% subplot(2,1,1)
% hold on
% plot(x, y, 'b', 'LineWidth', 1.5)
% patch(x, y, [0.5 0.5 0.5],'FaceColor','w');
% line([mu, mu],[floor(min(y)) max(y)],'Color','r','LineStyle','--')
% line([mu + sigma,mu + sigma],[0 mean(y)],'Color','k','LineStyle',':')
% line([mu + 2*sigma,mu + 2*sigma],[0 mean(y)],'Color','k','LineStyle',':')
% line([mu + 3*sigma,mu + 3*sigma],[0 mean(y)],'Color','k','LineStyle',':')
% text([mu + sigma,mu + sigma],[mean(y) mean(y)],' \sigma');
% text([mu + 2*sigma,mu + 2*sigma],[mean(y) mean(y)],' 2 \sigma');
% text([mu + 3*sigma,mu + 3*sigma],[mean(y) mean(y)],' 3 \sigma');
% %, num2str(mu,3)
% hold off
% 
% pd = makedist('Normal',mu,sigma);
% y2 = cdf(pd,x);
% 
% subplot(2,1,2)
% hold on
% plot(x,y2)
% hold off

%////////////// Extra Plotting / Images
t = linspace(0,2*pi);
xc = [const.d_small_actual/2 * cos(t); final_d_small_mu/2 * cos(t); const.d_big_actual/2 * cos(t); final_d_big_mu/2 * cos(t)];
yc = [const.d_small_actual/2 * sin(t); final_d_small_mu/2 * sin(t); const.d_big_actual/2 * sin(t); final_d_big_mu/2 * sin(t)];
figure('Name','Ring Plots','units','normalized','position',[.1 .1 .5 .8]);
hold on
plot(xc(1,:),yc(1,:),'r')
plot(xc(2,:),yc(2,:),'b')
plot(xc(3,:),yc(3,:),'r')
plot(xc(4,:),yc(4,:),'b')
legend({'Theory','Data'})
set(gca,'xtick',[],'ytick',[],'box','on')
hold off

% figure('Name','Plot Values')
% hold on
% rectangle('Position',[.5 0 1 final_d_small_mu])
% rectangle('Position',[2.5 0 1 final_d_big_mu])
% errorbar(1,final_d_small_mu,final_d_small_sigma)
% errorbar(3,final_d_big_mu,final_d_big_sigma)
% set(gca,'xlim',[0 4], 'ylim', [0 .22])
% hold off
figure('Name','BOX PLOTS BABY','units','normalized','position',[.1 .1 .7 .7]);
subplot(2,2,1)
boxplot(data_small_D')
subplot(2,2,2)
boxplot(data_small_d')
subplot(2,2,3)
boxplot(data_big_D)
subplot(2,2,4)
boxplot(data_big_d')

% zTest_weighted = std_small_D .* data_small_D;
% zTest_fit = polyfit(V,mean(zTest_weighted,2),1);
% zTest_val = polyval(zTest_fit,V);
% figure(99)
% plot(Vo,zTest_val)
%/////////////Clear Shit
clear t str xc yc
