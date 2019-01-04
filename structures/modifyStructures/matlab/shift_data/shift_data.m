clear all;

fid = fopen('data.relaxed');

tmp = textscan(fid, '%d%s', 1, 'CommentStyle', 'LAMMPS');
sys = tmp{1,1};

tmp = textscan(fid, '%s', 16);
tmp = textscan(fid, '%d%f', 2);
masses = tmp{1,2};

tmp = textscan(fid, '%s', 3);
tmp = textscan(fid, '%d%d%f%f%f%d%d%d', sys);

fclose(fid);

ids = double(tmp{1,1});
types = double(tmp{1,2});
x = tmp{1,3};
y = tmp{1,4};
z = tmp{1,5};
Pos = [ids, types, x, y, z];

[~,sort] = sort(Pos(:,3));
Pos = Pos(sort,:);

len = size(Pos,1)*(1/4);
ids = [];

for i = 1:len
    if Pos(i,2) == 2
        ids = [ids; i];
    end
end
        
block = Pos(ids,:);
Pos(ids,:) = [];

x_shift = (mean(Pos(end-128+1:end,3)) - mean(Pos(end-2*128+1:end-128,3)));
disp('x_shift =');
disp(x_shift);

x_min = min(Pos(:,3));
if x_min <= 0
	Pos(:,3) = (Pos(:,3) + abs(x_min));
else
	Pos(:,3) = (Pos(:,3) - abs(x_min));
end
disp('Pos min =');
disp(min(Pos(:,3)));
disp('Pos max =');
disp(max(Pos(:,3)));

block_min = mean(block(1:48,3));
%block_min = min(block(:,3));
if block_min <= 0
	block(:,3) = (block(:,3) + abs(block_min));
else
	block(:,3) = (block(:,3) - abs(block_min));
end
disp('Block min =');
disp(min(block(:,3)));

x_max = mean(Pos(end-128+1:end,3));
disp('Pos max =');
disp(x_max);

block(:,3) = (block(:,3) + x_max + x_shift);

Pos = [Pos; block];

%%
buffer = 0.009;

xlo = min(Pos(:,3))-buffer-20; xhi = max(Pos(:,3))+buffer+20;
ylo = min(Pos(:,4))-buffer; yhi = max(Pos(:,4))+buffer;
zlo = min(Pos(:,5))-buffer; zhi = max(Pos(:,5))+buffer;

%% Build the Header
fid = fopen(sprintf('data.%d.%s',sys,'shifted'),'w');
fprintf(fid,'LAMMPS data file\r\n\r\n%d atoms\r\n\r\n%d atom types\r\n\r\n'...
    ,sys,length(masses));
fprintf(fid,'%6.3f %6.3f xlo xhi\r\n%6.3f %6.3f ylo yhi\r\n%6.3f %6.3f zlo zhi\r\n'...
    ,xlo,xhi,ylo,yhi,zlo,zhi);
fprintf(fid,'\r\nMasses\r\n\r\n');
for i = 1:length(masses)
    fprintf(fid,'%d %6.4f\r\n',i,masses(i));
end
fprintf(fid,'\r\nAtoms\r\n\r\n');
%% Atoms
for i = 1:sys  
    fprintf(fid,'%d\t%d\t%6.3f\t%6.3f\t%6.3f\r\n',i,Pos(i,2),Pos(i,3),Pos(i,4),Pos(i,5));
end
fclose(fid);

%%
fid = fopen(sprintf('Si.Ge.%s.%d.xyz','shifted',sys),'w');
fprintf(fid,'%d\r\n',sys);
fprintf(fid,'%4.2f\t0\t0\t0\t%4.2f\t0\t0\t0\t%4.2f\r\n',max(Pos(:,3)),max(Pos(:,4))...
    ,max(Pos(:,5)));
for i = 1:sys
    if Pos(i,2)==1
        fprintf(fid,'Si\t%4.2f\t%4.2f\t%4.2f\r\n',Pos(i,3),Pos(i,4),Pos(i,5));
    elseif Pos(i,2)==2
        fprintf(fid,'Ge\t%4.2f\t%4.2f\t%4.2f\r\n',Pos(i,3),Pos(i,4),Pos(i,5));
    end
end
fclose(fid);






