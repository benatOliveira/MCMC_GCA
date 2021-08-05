
function [y,I]=mink(A,k)
[C,B]=sort(A);
I=B(1:k);
y=C(1:k);
end