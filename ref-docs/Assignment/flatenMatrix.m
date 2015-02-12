%maps repeating matrix values to unique values
function [A,cardinality]=flatenMatrix(A)
N=max(A);
vals=sort(A(:));
j=0;
prevVal=0;
for i=1:length(vals)
    if vals(i)>prevVal
        j=j+1;       
        prevVal=vals(i);
        I=find(A==vals(i));
        A(I)=j;
    end               
end 
cardinality=j; %the cardinality of the matrix or rather the size of the class
end