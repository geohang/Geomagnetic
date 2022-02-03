%this code is for Multiscale principal 
%component analysis of geomagnetic signal
% Hang Chen Feb 2, 2022

clear

load('CNB2000-2004.mat')

% example fro Z direction signal
 Z_mag = detrend(BZ1);
 x = Z_mag;
level = 5;
wname = 'sym4';
npc = 'heur';
[x_sim,qual,npcA] = wmspca(Z_mag,level,wname,npc);


figure 
subplot(1,2,1)
plot(x)
title('Original signal')
subplot(1,2,2)
plot(x_sim)
title('PCA signal')

npcA(1:3) = zeros(1,3);
[x_sim,qual,npcB] = wmspca(x,level,wname,npcA);
