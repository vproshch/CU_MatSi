function Poso = SurfaceWelds2(Pos,a)
Pos = RemoveRepeats(Pos);
Pos = Pos - [min(Pos(1,:));min(Pos(2,:));0];
posx = Pos(1,:);
xmax = max(posx);
nx = round(xmax/a);
shift = a/10;

for i = 0:.5:nx
    yzl = Pos(:,find(posx == a*i));
    yzu = Pos(:,find(posx == a*(i+.25)));
    zmax = max(yzu(3,:));
    zmin = min(yzl(3,:));
    yz_upper = yzu(:,find(yzu(3,:)==zmax))
    yz_lower = yzl(:,find(yzl(3,:)==zmin))
    for j = 1:length(yz_upper)
        if rem(yz_upper(1,j)-yz_upper(2,j)-a/2,2*a)
            s = -shift;
        else 
            s = shift;
        end
        yz_upper(1,j) = yz_upper(1,j) + s;
        yz_upper(2,j) = yz_upper(2,j) - s;
    end
    for j = 1:length(yz_lower)
        if rem(yz_lower(1,j)-yz_lower(2,j),2*a)
            s = -shift;
        else 
            s = shift;
        end
        yz_lower(1,j) = yz_lower(1,j) + s;
        yz_lower(2,j) = yz_lower(2,j) - s;
    end
    yzu = yzu(:,find(yzu(3,:)~=zmax));
    yzl = yzl(:,find(yzl(3,:)~=zmin));
    yzu = [yzu yz_upper]; yzl = [yzl yz_lower]
    Pos(:,find(posx == a*i)) = yzl;
    Pos(:,find(posx == a*(i+.25))) = yzu;
end
Poso = Pos;
end
        
    
            