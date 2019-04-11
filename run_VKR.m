close all
clear all

load data.mat

% prepare input time series

input{1}=Precip;
input{2}=PET;

% set parameters

param(1)=4.7; % a
param(2)=45; % K_epi
param(3)=1200; % S_soil
param(4)=260; % S_epi

% run model

[sim]=VarKarst_R(input,param);

R=sim(:,1);
AET=sim(:,2);

figure;subplot(2,1,1)
plot(date(:,4),[Precip R]); datetick('x',10);legend('Prec obs [mm/d]','R sim [mm/d]');
subplot(2,1,2);plot(date(:,4),AET_obs,'x'); 
hold all;plot(date(:,4),AET); 
datetick('x',10);legend('AET obs [mm/d]','AET sim [mm/d]');