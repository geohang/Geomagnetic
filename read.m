%%
clear;
num1=xlsread('C:\Users\lenovo\Desktop\CNB\2000.xlsx');
num2=xlsread('C:\Users\lenovo\Desktop\CNB\2001.xlsx');
num3=xlsread('C:\Users\lenovo\Desktop\CNB\2002.xlsx');
num4=xlsread('C:\Users\lenovo\Desktop\CNB\2003.xlsx');
num5=xlsread('C:\Users\lenovo\Desktop\CNB\2004.xlsx');
num6=xlsread('C:\Users\lenovo\Desktop\CNB\2005.xlsx');
num7=xlsread('C:\Users\lenovo\Desktop\CNB\2006.xlsx');
num8=xlsread('C:\Users\lenovo\Desktop\CNB\2007.xlsx');
num9=xlsread('C:\Users\lenovo\Desktop\CNB\2008.xlsx');
num10=xlsread('C:\Users\lenovo\Desktop\CNB\2009.xlsx');
num11=xlsread('C:\Users\lenovo\Desktop\CNB\2010.xlsx');
num12=xlsread('C:\Users\lenovo\Desktop\CNB\2011.xlsx');
num13=xlsread('C:\Users\lenovo\Desktop\CNB\2012.xlsx');
num14=xlsread('C:\Users\lenovo\Desktop\CNB\2013.xlsx');
num15=xlsread('C:\Users\lenovo\Desktop\CNB\2014.xlsx');
num16=xlsread('C:\Users\lenovo\Desktop\CNB\2015.xlsx');
num17=xlsread('C:\Users\lenovo\Desktop\CNB\2016.xlsx');
num18=xlsread('C:\Users\lenovo\Desktop\CNB\2017.xlsx');
num19=xlsread('C:\Users\lenovo\Desktop\CNB\2018.xlsx');
num20=xlsread('C:\Users\lenovo\Desktop\CNB\2019.xlsx');
num21=xlsread('C:\Users\lenovo\Desktop\CNB\2020.xlsx');
num22=xlsread('C:\Users\lenovo\Desktop\CNB\2021-2022.xlsx');

%%
% ALL=[num1;num2;num3;num4;num5;num6;num7;num8;num9;num10;num11;num12;num13;num14;num15;num16;num17;num18;num19;num20;num21;num22;];
% H_N=ALL(:,1);
% BZ1=ALL(:,3);
% save('CNB2000-2022.1.20.mat','H_N','BZ1');
ALL=[num6;num7;num8;num9;num10;];
H_N=ALL(:,1);
BZ1=ALL(:,3);
save('CNB2005-2010.mat','H_N','BZ1');