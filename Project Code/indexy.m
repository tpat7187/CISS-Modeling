function binarytodecimal = indexy(A)
s = 0;
for j = 1:length(A)
        s = s + (2^(length(A) - j))*A(j);
    end
    binarytodecimal = s;
end
