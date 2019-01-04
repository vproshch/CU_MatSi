%{
This functions takes a matrix of Si atomic positions as well as a string name
and writes a .xyz file for VMD visualization.

    INPUTS
        Pos: Is a 3xN matrix where Pos(1:3,i) is the x y and z coordinates
            of the ith atom
        name: Is a string that determines the name of the output file
    OUTPUTS
        This function writes a file Si.name.N.xyz that can be called for
        VMD visuallization
%}
function writeXYZ(Pos,name)
sys = length(Pos);
fid = fopen(sprintf('Si.%s.%d.xyz',name,sys),'w');
fprintf(fid,'%d\r\n',sys);
fprintf(fid,'%4.2f\t0\t0\t0\t%4.2f\t0\t0\t0\t%4.2f\r\n',max(Pos(1,:)),max(Pos(2,:))...
    ,max(Pos(3,:)));
for i = 1:sys
    fprintf(fid,'Si\t%4.2f\t%4.2f\t%4.2f\r\n',Pos(1,i),Pos(2,i),Pos(3,i));
end
fclose(fid);
end