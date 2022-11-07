function id = getEqId(i,j)
global flag gEq;

id = -flag(i,j);
if (id < 0)
    id = gEq(flag(i,j));
end
