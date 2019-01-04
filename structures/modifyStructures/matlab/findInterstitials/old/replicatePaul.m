close all; clear all; clc;

%%% THIS SCRIPT TAKES A MATLAB WORKSPACE GENERATED FROM THE "makeBulkPaul"
%%% ENTER THE NAME OF THIS WORKSPACE FILE DOWN BELOW (THE "load" LINE)
%%% THIS WORKSPACE MUST CONTAIN A "Pos" VECTOR CONTAINING ALL THE BULK
%%% ATOMS, AS WELL AS 1 INTERSTITIAL SITE AT END OF THE VECTOR.

%%% THIS SCRIPT REPLICATES THE INTERSTITIAL SITE THROUGHOUT THE STRUCUTRE
%%% BY A LATTICE CONSTANT IN ALL DIRECTIONS, THEN RANDOMLY POPULATES N
%%% NUMBER OF THESE SITES BASED ON HOW MANY YOU WANT. AT END, WRITES AN XYZ
%%% FILE AND A TRJ FILE FOR LAMMPS AND VMD VIEWING

%%% NEXT STEPS: RELAX THIS STRUCTURE IN LAMMPS, THEN AFTER RELAXATION YOU
%%% CAN PLOT THE STRUCUTRE AND ANALYZE NEAREST NEIGHBORS IN THE "AFTER
%%% RELAX PLOT" SCRIPT. YOU CAN USE "readLAMMPS" SCRIPT TO GENERATE A POS
%%% VECTOR FROM A LAMMPS DATA FILE


%%INPUTS:
numAdd = 64; %number of interstitial atoms to add (HERE IS WHERE YOU SET THIS)
save_file_name = '8x8x8_64int'; %in here, put what you want the output files to be called

%-----------------------------------------------------
atom_size = 6; %just used for plotting... ignore this

%load in a matlab file with just 1 interstitial already added (which is in the
%last elemnet of Pos matrix)
load('workspace.mat');         %load('point_1_mesh.mat')
%plot the lattice atoms and the one interstitial site
figure(1)
plot3(Pos(1,1:end-1),Pos(2,1:end-1),Pos(3,1:end-1),'o','MarkerEdgeColor','k','MarkerFaceColor',...
    'g','MarkerSize',atom_size)
hold on;
plot3(Pos(1,end),Pos(2,end),Pos(3,end),'*','MarkerEdgeColor','r','MarkerFaceColor',...
    'r','MarkerSize',atom_size/2)

figure(2) %figure for first interstitial atom distances
scatter(1:20,intDists(p,1:20));
xlabel('atom #');
ylabel('distance (A)');
title('interstitial nearest neighbor distances');

%%%% HERE we will replicate the original interstitial site to fill the
%%%% whole region

%replicating to find possible spots
potential_sites = Pos(:,end);
original = potential_sites;

    %x direction (producing a 'line' of atoms
    numRepeat = 30; %number of times to replicate in each direction. can be increased if needed
    for i = 1:numRepeat
        
        add_forward = original + [i*a 0 0]';
        add_back = original - [i*a 0 0]';
        
        if (add_forward(1) < xhi) %if trying to add outside bounds, automatically end the loop
            potential_sites(:,end+1)  = add_forward;
        end
        
        if (add_back(1) > xlo)
            potential_sites(:,end+1)  = add_back;
        end
       
    end
    
%y direction (just duplicate this 'line' of atoms)
    %at this point, potential_sites contains the line of #columns = #
    %atoms
    atomLine = potential_sites;

for i = 1:numRepeat
    line_forward = (atomLine + [0 i*a 0]'.*ones(size(atomLine))); %shift this line and add to list of sites
    line_back = (atomLine - [0 i*a 0]'.*ones(size(atomLine)));

    %check if lines are out of bounds
    if max(line_forward(2,:)) < yhi
        potential_sites = [potential_sites line_forward];
    end
     if min(line_back(2,:)) > ylo
        potential_sites = [potential_sites line_back];
    end

end
        
%z direction:
    %at this point, potential_sites contains a "layer" of atoms
    atomLayer = potential_sites;

for i = 1:numRepeat
    %replicate layer
    layer_forward = (atomLayer + [0 0 i*a]'.*ones(size(atomLayer))); %shift this line and add to list of sites
    layer_back = (atomLayer - [0 0 i*a]'.*ones(size(atomLayer)));
    
    %check if layers exceed bounds
    if max(layer_forward(3,:)) < zhi
        potential_sites = [potential_sites layer_forward];
    end
     if min(layer_back(3,:)) > zlo
        potential_sites = [potential_sites layer_back];
    end

end
        
%plotting all possible sites
figure(1);
hold on;
plot3(potential_sites(1,:),potential_sites(2,:),potential_sites(3,:),'*','MarkerEdgeColor','r','MarkerFaceColor',...
    'r','MarkerSize',atom_size/2)

%%%%%%%%%%%%%
%at this point, we have a list of all potential_sites replicated from the
%original interstitial. now we pick random spots to add an atom

%total number of possibe interstitial sites
numSites = length(potential_sites);


interstitials = [];
indexList = []; %track index of sites already filled
for i = 1:numAdd
    
    addedFlag = 0;
    
    while addedFlag == 0 %try to pick a random spot to add interstitial. keep trying until unique one found
        index = randi(numSites);
        
        if isempty(find(indexList==index,1))
            addAtom = potential_sites(:,index);
            addedFlag = 1;
            indexList = [indexList index]; %add the index to the list, so we dont pick same spot twice
        end
    
    end
    
    interstitials = [interstitials addAtom];
    
end

%plot the filled interstitial sites
hold on;
plot3(interstitials(1,:),interstitials(2,:),interstitials(3,:),'o','MarkerEdgeColor','y','MarkerFaceColor',...
    'y','MarkerSize',atom_size)
        
%preparing Pos vector to write to data file
finalPositions = [Pos(:,1:end-1) interstitials];
%shift everything by -0.2
finalPositions = finalPositions - 0.2;

finalPositions(4,:) = ones(1,length(finalPositions)); %this 4th row specifies atom group number for lammps (all 1's if all silicon)



%write lammps data file
writeLAMMPS(finalPositions,save_file_name,28.0855,0,0,0);

%also write xyz file for VMD viewing
writeXYZ(finalPositions(1:3,:),save_file_name);

Pos = [Pos(:,1:end-1) interstitials];
inters = interstitials;
numInt = numAdd;
%analyzing nearest distances
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
        
    if p < 4 %only plot for this range 
        figure(1+p)
        scatter(1:20,intDists(p,1:20));
        xlabel('atom #');
        ylabel('distance (A)');
        title('interstitial nearest neighbor distances');
    end
    
    %getting number of nearest neighbors
    jumps(p,:) = diff(intDists(p,:));
    index = find (jumps(p,:) > 0.3,1);
    numNear(p) = index;
    
end

figure
hist(numNear)



function [distance]= getDistance(r,atomPos)
    distance =( ( r(1) - atomPos(1) )^2 + ( r(2) - atomPos(2) )^2 + ( r(3) - atomPos(3) )^2 )^(1/2);
end
