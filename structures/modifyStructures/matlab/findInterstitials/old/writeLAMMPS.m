%{
This functions takes a matrix of Si atomic positions as well as a string name
and other information and writes a LAMMPS data file.

    INPUTS
        Pos: Is a 4xN matrix where Pos(1:3,i) is the x y and z coordinates
            of the ith atom and Pos(4,i) is a whole number that determines
            the atom type.  If all of the atoms are the same, then this row
            should be all ones.  In this case, the output of the functions 
            makeGrooves, readXYZ or readLAMMPS needs to be modified via the
            following line of code to be compatible with this function.
                >> Pos = [Pos;ones(1,length(Pos))];
        name: Is a string that determines the name of the output file
        masses: A vector of the molecular weights of the atoms where
                masses(n) cooresponds to the molecular weight of the of any
                atom that has Pos(4,i) == n.  If there is one atom type,
                masses is a scalar.
        xbuffer: Vacuum buffer in the x direction in angstroms on either
                side of the structure
        ybuffer: Vacuum buffer in the y direction in angstroms on either
                side of the structure
        zbuffer: Vacuum buffer in the z direction in angstroms on either
                side of the structure
    OUTPUTS
        This function writes a file data.bulk.Si.N.name to read in as a data
        file for LAMMPS
%}

function writeLAMMPS(Pos,name,masses,xbuffer,ybuffer,zbuffer)
xlo = min(Pos(1,:))-xbuffer; xhi = max(Pos(1,:))+xbuffer;
ylo = min(Pos(2,:))-ybuffer; yhi = max(Pos(2,:))+ybuffer;
zlo = min(Pos(3,:))-zbuffer; zhi = max(Pos(3,:))+zbuffer;
sys = length(Pos);
%% Build the Header
fid = fopen(sprintf('data.bulk.Si.%d.%s',sys,name),'w');
fprintf(fid,'LAMMPS data file for Si structure\r\n\r\n%d atoms\r\n\r\n%d atom types\r\n\r\n'...
    ,sys,length(masses))
fprintf(fid,'%6.3f %6.3f xlo xhi\r\n%6.3f %6.3f ylo yhi\r\n%6.3f %6.3f zlo zhi\r\n'...
    ,xlo,xhi,ylo,yhi,zlo,zhi)
fprintf(fid,'\r\nMasses\r\n\r\n')
for i = 1:length(masses)
    fprintf(fid,'%d %6.4f\r\n',i,masses(i))
end
fprintf(fid,'\r\nAtoms\r\n\r\n')
%% Atoms
for i = 1:sys  
    fprintf(fid,'%d\t%d\t%6.3f\t%6.3f\t%6.3f\r\n',i,Pos(4,i),Pos(1,i),Pos(2,i),Pos(3,i));
end
fclose(fid);
end