function [Pos_Vector] = makeBulkPaul(nx,ny,nz,nzi,nzj)
%%% THIS FUNCTION CREATES A BULK SILICON STRUCTURE OF NX by NY by NZ UNIT CELLS. 
%%% (Make NZI and NZJ = to 0 if no grooves are desired)

%%% THEN, THIS FUNCTION LOCATES AN INTERSTITIAL SITE USING GEOMETRIC
%%% MINIMIZATION ALGORITHM

%%%OUTPUTS: THE FUNCTION WILL OUTPUT A POSITION VECTOR CONTAINING XYZ
%%%COORDINATES OF ALL ATOMS, AND THE LAST ELEMENT WILL BE THE POTENTIAL
%%%INTERSTITIAL SITE

%%%NEXT STEP: SET A BREAKPOINT AT END OF THIS FUNCTION, RUN IT TO GENERATE
%%%YOUR STRUCTURE, THEN SAVE THE WORKSPACE AS A .MAT FILE. THEN, GO TO THE
%%%"replicatePaul" SCRIPT TO PROCEED.



tic
close all;
sys = nx*ny*nz*8 + length(2:2:nx)*ny*8; 
a = 5.431;
fid = fopen(sprintf('data.bulk.Si.%d.grooves.txt',sys),'w');

%Heading (these should not be fixed probably?)
% fprintf(fid,'LAMMPS Data file for Si with grooves\r\n\r\n%d atoms\r\n\r\n1 atom type\r\n\r\n',sys);
% xlo = -.2; xhi = nx*a;
% ylo = -.2; yhi = ny*a;
% zlo = -.2; zhi = (nz+1)*a;
% fprintf(fid,'%4.4f %4.4f xlo xhi\r\n%4.4f %4.4f ylo yhi\r\n%4.4f %4.4f zlo zhi\r\n\r\n',...
%     xlo,xhi,ylo,yhi,zlo,zhi);
% fprintf(fid,'Masses\r\n\r\n1 28.0855\r\n\r\nAtoms\r\n\r\n');

%The atoms
xi = 0;
yi = 0;
zi = 0;
xst = 0;
yst = 0;
zst = 0;
st = -nzj/2;
n = 1;
Posx = [];
for i = 1:8
    for j = 1:nx
        for k = 1:ny
            for l = 1:(nz+nzj/2+st)
                fprintf(fid,'%-d\t1\t%-4.4f\t%-4.4f\t%-4.4f\r\n',...
                    n,xi+xst,yi+yst,zi+zst);
                Pos(:,n) = [xi + xst;yi + yst;zi + zst];
                n = n+1;
                zi = zi + a;
            end
            zi = 0;
            if round(k/2)==k/2
                if j > ny/2
                    zst = zst + a/2;
                else
                    zst = zst - a/2;
                end
            else
                if j > ny/2
                    zst = zst - a/2;
                else
                    zst = zst + a/2;
                end
            end
            yi = yi + a/2;
        end
        if round(j/2)==j/2
            yi = 0;
            yst = yst + a/4;
            xi = xi + a/4;
            zst = zst + a/4;
            if round(j/4) == j/4
                yst = 0;
            end
        end
    end
    if rem(i,nzi) == 0
        st = st*-1;
    end
    zst = 0;
end
figure
plot3(Pos(1,:),Pos(2,:),Pos(3,:),'o','MarkerEdgeColor','k','MarkerFaceColor',...
    'g','MarkerSize',10)
%axis([0 a 0 a 0 a]


%%%%%%%%%%%%%%%%CODE ADDED BY PAUL!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USING METHOD:
% https://www.sciencedirect.com/science/article/pii/S1359646207008871 

%finding bounds of domain
xlo = min(Pos(1,:));
xhi = max(Pos(1,:));
ylo = min(Pos(2,:));
yhi = max(Pos(2,:));
zlo = min(Pos(3,:));
zhi = max(Pos(3,:));

%discretizing the domain so we can check the potential at various
%positions, check only in interior    
dx = 0.1; dy = 0.1; dz = 0.1;
weight = 6; 
buffer = 18 %6.5; %buffer to only look in interior. 
%%% This is so there are fewer locations to search, and so interstitial are
%%% not always placed at very edge of the structure.


xVec = xlo+buffer:dx:xhi-buffer; %vector of all the coordinates to test
yVec = ylo+buffer:dy:yhi-buffer;
zVec = zlo+buffer:dz:zhi-buffer;



numInter = 1; %just do 1 iteration

for p = 1:numInter %for each interstitial site to find
    for x = 1:length(xVec) %loop through all xyz coords
       for y = 1:length(yVec)
          for z = 1:length(zVec)
              %executing the formula from the paper
              geoSum = 0;
              for i = 1:(n-1) %for every atom
                  %calculate distance to atom
                  dist = getDistance([xVec(x) yVec(y) zVec(z)],Pos(:,i) );
                  geoSum = geoSum + dist^(-weight);

              end
                S(x,y,z) = geoSum^(-1/weight); %the function we want to maximize

          end
       end
    end

     [~,maxSpot] = max(S(:)); %get 1-D array index that maximizes function S
     [xInd,yInd,zInd]=  ind2sub(size(S),maxSpot); %get this index as 3D indices
     figure(1)
     hold on;
     plot3(xVec(xInd),yVec(yInd),zVec(zInd),'o','MarkerEdgeColor','k','MarkerFaceColor',...
        'r','MarkerSize',10)

    %after we have found/visualized the interstitial site, add site to Pos vector 
    Pos(:,end+1) = [xVec(xInd); yVec(yInd); zVec(zInd)];  
    intSpot = Pos(:,end);
    n=n+1; %increment number of atoms
    %now we can repeat the process
    
    %sanity check, distance to nearest atoms
     for i = 1:(n-1 -1) %for every atom excluding itself
          %calculate distance to atom
          intDists(p,i) = getDistance(intSpot,Pos(:,i) );
          

     end
     intDists(p,:) = sort(intDists(p,:));
     avgDist = mean(intDists(p,1:8));
     figure(1+p)
     scatter(1:20,intDists(p,1:20));
     xlabel('atom #');
     ylabel('distance (A)');
     title('interstitial nearest neighbor distances');
end



function [distance]= getDistance(r,atomPos)
    distance =( ( r(1) - atomPos(1) )^2 + ( r(2) - atomPos(2) )^2 + ( r(3) - atomPos(3) )^2 )^(1/2);
end


toc

Pos_Vector = Pos;

save('workspace');

end

