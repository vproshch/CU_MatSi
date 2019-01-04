%{
This function reads a LAMMPS data file and converts it to a matrix of
postions

    INPUTS
        filename: string of filename to read

    OUTPUTS
        Pos: 3xN matrix where Pos(1:3,i) is the x y and z coordinates of the
            ith atom 
%}
function Pos = readLAMMPS(filename)
fid = fopen(filename);
%{
For a future revision
    trash = fgetl(fid);trash = fgetl(fid);
    numAtoms = fgetl(fid); numAtoms(end-5:end) = '';numAtoms = str2double(numAtoms);
    trash = fgetl(fid);trash = fgetl(fid);trash = fgetl(fid);
    xlim = fgetl(fid);xlim(end-7:end)='';C = strsplit(xlim,' ');xlim = str2double(C{2});
    ylim = fgetl(fid);ylim(end-7:end)='';C = strsplit(ylim,' ');ylim = str2double(C{2});
    zlim = fgetl(fid);zlim(end-7:end)='';C = strsplit(zlim,' ');zlim = str2double(C{2});
%}
for i = 1:9
    trash = fgetl(fid);
end
Pos = []; last = 1; i = 1;
while last
    pi = fgetl(fid)
    if pi == -1
        last = 0;
    elseif strcmp(pi,'')
        last = 0;
    else
        pi = regexprep(pi,'+','');
        pos = strsplit(pi,'\s*','DelimiterType','RegularExpression')
        Pos(1,i)= str2double(pos{3});
        Pos(2,i)= str2double(pos{4});
        Pos(3,i)= str2double(pos{5});
        i = i + 1
    end
end
fclose all
end