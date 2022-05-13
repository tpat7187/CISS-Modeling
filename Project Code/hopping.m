function [a,H] = hopping(a,b,c,spin1,spin2,cat_index,H,constant)
v_sight = 1:2:20;
v_sight2 = 2:2:20;

%calculates the position of the electrons on-site and in relation to the
%input state
sight1 = a(v_sight(b):v_sight(b)+1);
sight2 = a(v_sight(c):v_sight(c)+1);


if spin1 == "up" 
    position1 = sight1(1);
    var = v_sight(b);
    val = 1;
else 
    position1 = sight1(2);
    var = v_sight2(b);
    val = 2;
end
if spin2 == "up" 
    position2 = sight2(1);
    var2 = v_sight(c);
    val2 = 3;
else
    position2 = sight2(2);
    var2 = v_sight2(c);
    val2 = 4;
end



if position1 == 1 && position2 == 0
    index_pos1 = indexy(a);
    pos1 = find(cat_index == index_pos1);
    a(var) = 0;
    a(var2) = 1;
    if b < c
         sign = (-1)^(nnz(a(var:var2))-1);
    end
    if b > c
        sign = (-1)^(nnz(a(var2:var))-1);
    end
    index_pos2 = indexy(a);
    pos2 = find(cat_index == index_pos2);
    H(pos1,pos2) = sign*constant;
else
    sign = 0;
    pos1 = 0;
    pos2 = 0;
end
end

