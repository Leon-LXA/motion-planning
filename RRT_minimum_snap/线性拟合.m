clc;clear;close all;
x = [0.2 0.32 0.91 1.2];
y = [5.198 5.492 6.884 7.413];
plot(x,y,'color',[0,1.0,0],'LineWidth',2);
hold on;
y1 = [5:0.001:8];
c = polyfit(y,x,1);
x1 = y1*c(1)+c(2);
plot(x1,y1,'*');