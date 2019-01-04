function [Pos_Vector] = makeBulkPaul(nx,ny,nz)
%some stuff used in the code for grooves, but we don't care aboout.
ni = nx;
nI = nx;
nj = 0;

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

%quentins code
% INPUTS 
%         nx: The number of unit cells in the x direction. nx must be related
%             to ni and nI by nx = a*ni + b*nI where a and b are whole 
%             numbers
%         ny: The number of unit cells in the y direction nz: The number of 
%             unit cells in the z direction 
%         ni: The number of unit cells for the width of the grooves. 
%         nI: The number of unit cells between grooves 
%         nj: The number of unit cells for the height of the grooves both 
%             above and below the main block of atoms
%         
% OUTPUTS Pos: matrix of atom positions.  Pos(:,i) is a 3 x 1 vector giving
% the x y and z position of the ith atom.
% 
% If a block of atoms is desired without any grooves then set ni = nx
% %}

b = 4;%Used to build the structure to aviod rounding error
a = 5.431;%Lattice Constant
%{
Unit = [0 0  0  0 0 a/2 a/2 a/2 a/2 a a  a  a a a/4 3*a/4   a/4 3*a/4;...
        0 0 a/2 a a  0  a/2 a/2  a  0 0 a/2 a a a/4 3*a/4 3*a/4   a/4;...
        0 a a/2 0 a a/2  0   a  a/2 0 a a/2 0 a a/4   a/4 3*a/4 3*a/4;];
 %}
%The unit cell
Unit = [0  0  b/2 b/2 b/4 3*b/4   b/4 3*b/4;...
        0 b/2  0  b/2 b/4 3*b/4 3*b/4   b/4;...
        0 b/2 b/2  0  b/4  b/4  3*b/4 3*b/4;];
%plot3(Unit(1,:),Unit(2,:),Unit(3,:),'.','MarkerSize',30)    
%Replicate in the x direction
Sm = [];
Bg = [];
for i = 0:nI-1
    Sm = [Sm Unit+b.*[i;0;0]];
end
for i = 0:ni-1
    Bg = [Bg Unit+b.*[i;0;0]];
end
Sm = RemoveRepeats(Sm);
Bg = RemoveRepeats(Bg);
%Replicate in the y direction
PosjS = Sm;
PosjB = Bg;
for j = 0:ny-1
    Sm = [Sm PosjS+b.*[0;j;0]];
    Bg = [Bg PosjB+b.*[0;j;0]];
end
Sm = RemoveRepeats(Sm);
Bg = RemoveRepeats(Bg);
%replicate in the z direction for the main block height
Posk = Sm;
for k = 0:nz-1
    Sm = [Sm Posk+b.*[0;0;k]];
end
Sm = RemoveRepeats(Sm);
%replicate in the z direction for the grooves height
Posk = Bg;
for k = -nj:nz+nj-1
    Bg = [Bg Posk+b.*[0;0;k]];
end
Bg = RemoveRepeats(Bg);

Pos = [];
num_steps = 2*floor(nx/(nI+ni))+ ceil(rem(nx,ni+nI)/nI);
for i = 0:num_steps-1
    if i/2 == round(i/2)
        Pos = [Pos Sm+b*(ni+nI)*(i/2).*[1;0;0]];
    else
        Pos = [Pos Bg+b*(nI*(i+1)/2+ni*(i-1)/2).*[1;0;0]];
    end
end
Pos = Pos + abs(min(Pos(3,:))).*[0;0;1];%Shift the structure so all z positions are positive
Pos = RemoveRepeats(Pos);

%plot3(Pos(1,:),Pos(2,:),Pos(3,:),'.')
%Pos = SurfaceWelds2(Pos,b); % Add diamurization (not sure if that is spelled correctly)
Pos = Pos.*(a/b);%Convert to lattice constant lengths
plot3(Pos(1,:),Pos(2,:),Pos(3,:),'o','MarkerEdgeColor','k','MarkerFaceColor',...
    'g','MarkerSize',10)


n = length(Pos) + 1;
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

axis([0 1.2*xhi 0 1.2*xhi 0 1.2*xhi]);

%discretizing the domain so we can check the potential at various
%positions, check only in interior    
dx = 0.1; dy = 0.1; dz = 0.1;
weight = 6; 
buffer = 7%18 %6.5; %buffer to only look in interior. 
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
save('workspace.mat')
end

