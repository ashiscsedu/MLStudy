function e=classification_error(s,I)
%attempt to calculate classification error in percentage.
%should it give the same answer as the Hopkins155 code?

k=max(s); %determine the number of objects from the ground truth vector

l=length(I);
Sc=zeros(1,k);
P=perms(1:k); %Permutation matrix
%loop over all the k permutations
for i=1:size(P,1)
    IP=zeros(1,l);
    for j=1:l
        for n=1:k
           if I(j)==n
               IP(j)=P(i,n);
           end
        end
    
        
    end
    
    %test current permutation with ground truth s.
    
    Sc(i)=sum(s==IP');
end

e=(1-max(Sc)/length(s))*100;

end