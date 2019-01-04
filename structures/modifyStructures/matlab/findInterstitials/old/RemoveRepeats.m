function Out = RemoveRepeats(In)
n = 1;
In = round(In.*10000)./10000;
while n < length(In)
    vec = round(In(:,n));
    rm = [];
    if n == 136
        j=0
    end
    for i = n+1:length(In)
        if round(In(1,i)) == vec(1) && round(In(2,i)) == vec(2) && round(In(3,i)) == vec(3)
            rm = [rm i];
        end
    end
    In(:,rm) = [];
    n = n + 1;
end
Out = In;
end