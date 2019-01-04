%{
This function generates a .xyz file for vmd visualization and an input data
file for LAMMPS of an fcc Silcon structure with grooves. nx, ny and nz
define the dimensions of the main block while ni and nj define the width
and height of the grooves.

    INPUTS 
        nx: The number of unit cells in the x direction. nx must be related
            to ni and nI by nx = a*ni + b*nI where a and b are whole 
            numbers
        ny: The number of unit cells in the y direction 
        nz: The number of unit cells in the z direction 
        ni: The number of unit cells for the width of the grooves. 
        nI: The number of unit cells between grooves 
        nj: The number of unit cells for the height of the grooves both 
            above and below the main block of atoms
        write: logical that evaluates as true if the .xyz file and Lammps
            input file are desired
        *********************
        TO MAKE BULK MEMEBRANE: nx, ny, and nz are desired dimensions. Set ni
        and nj = 1 to keep from causing problems with loops, but set nI =
        nx. i.e. distance between structures = lenght of membrane, so no
        structures are made.
        *********************
OUTPUTS Pos: matrix of atom positions.  Pos(:,i) is a 3 x 1 vector giving
the x y and z position of the ith atom.

If a block of atoms is desired without any grooves then set nI = nx, ni and
nj not equal to 0.
%}
function Pos = makeGrooves(nx,ny,nz,ni,nI,nj,write)
b = 4;%Used to build the structure to aviod rounding error
a = 3.570;%Lattice Constant #5.431 Si
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

plot3(Pos(1,:),Pos(2,:),Pos(3,:),'.')
Pos = SurfaceWelds2(Pos,b); % Add diamurization (not sure if that is spelled correctly)
plot3(Pos(1,:),Pos(2,:),Pos(3,:),'.')
Pos = Pos.*(a/b)%Convert to lattice constant lengths
%% Write xyz and lammps file
if write
    str = input('Name of Structure:\n','s');
    writeXYZ(Pos,str);
    Pos = [Pos;ones(1,length(Pos))];
    xbuffer = input('Vacuum in the x direction (A):\n');
    ybuffer = input('\nVacuum in the y direction (A):\n');
    zbuffer = input('\nVacuum in the z direction (A):\n');
    writeLAMMPS(Pos,str,28.0885,xbuffer,ybuffer,zbuffer)
end
end
        