%%% THIS SCRIPT TAKES A MATLAB WORKSPACE CONTAINING JUST THE POS VECTOR OF
%%% ATOMS (ORIGINALS ALL FIRST, AND THE LAST ELEMENTS ARE THE
%%% INTERSTITIALS)

%%% THIS SCRIPT PLOTS THE RELAXED STRUCTURE, AND GIVES INFORMATION ABOUT
%%% THE NEAREST NEIGHBORS


%plotting after relax
clear all; close all; clc;
atomsize = 6;
load('after_relax_positions.mat');
Pos = ans;

sys = 4096; %make sure this matches the number of atoms excluding interstitials
numInt = 64; %make sure this matches the number of interstitials

figure(1)
plot3(Pos(1,1:sys),Pos(2,1:sys),Pos(3,1:sys),'o','MarkerEdgeColor','k','MarkerFaceColor',...
    'g','MarkerSize',atomsize)
hold on;
plot3(Pos(1,(sys+1):end),Pos(2,(sys+1):end),Pos(3,(sys+1):end),'o','MarkerEdgeColor','k','MarkerFaceColor',...
    'r','MarkerSize',atomsize)

%analyzing closest distances
originals = Pos(:,1:sys);
inters = Pos(:,sys+1:sys+numInt);

Pos = [originals inters];
for p = 1:numInt
intSpot = inters(:,p);

    for i = 1:(length(Pos)) %for every atom excluding itself
      %calculate distance to atom
      if intSpot ~= Pos(:,i)
        intDists(p,i) = getDistance(intSpot,Pos(:,i));
      else
        intDists(p,i)= nan;
      end
      
    end
    
    
        intDists(p,:) = sort(intDists(p,:));
        
    if p < 20 %CHANGE THIS TO PLOT MORE OR FEWER
        figure(1+p)
        scatter(1:20,intDists(p,1:20));
        xlabel('atom #');
        ylabel('distance (A)');
        title('interstitial nearest neighbor distances');
    end
    
    %getting number of nearest neighbors
    jumps(p,:) = diff(intDists(p,:));
    index = find (jumps(p,:) > 0.2,1);
    numNear(p) = index;
    
end

figure
histogram(numNear)
xlabel('number of nearest neighbors');
ylabel('count');
title('nearest neighbors of the 64 interstitial atoms');

function [distance]= getDistance(r,atomPos)
    distance =( ( r(1) - atomPos(1) )^2 + ( r(2) - atomPos(2) )^2 + ( r(3) - atomPos(3) )^2 )^(1/2);
end