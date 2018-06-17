clc;clear;close all;

Data1 = xlsread('0102selected.xlsx');
Data2 = xlsread('2807selected.xlsx');

epsilong0 = 1;
[von_result1,lamda1]=vonfun(log(Data1(:,2)),epsilong0,size(Data1,1));
[von_result2,lamda2]=vonfun(log(Data2(:,2)),epsilong0,size(Data2,1));

fig1 = figure(1);
plot(log(Data1(:,1)),log(Data1(:,2)),log(Data1(:,1)),von_result1)
title('e');
legend('1','2');
fig2 = figure(2);
plot(log(Data2(:,1)),log(Data2(:,2)),log(Data2(:,1)),von_result2)
title('s');
legend('1','2');