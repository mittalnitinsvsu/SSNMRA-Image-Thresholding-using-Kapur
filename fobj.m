%% fitness function for Kapur Entropy
function fnew=fobj(u,nd,probR)
PI0=zeros(1,size(probR,2));
for i=1:size(probR,2)
PI0(:,i) = probR(:,i); 
    ind = PI0(i) == 0;
    ind = ind .* eps;
    PI0(:,i) = PI0(:,i) + ind;
end
probR=PI0;
j=1;
w1=sum(probR(1:u(j,1)));
n1=(probR/w1).*log(probR/w1);
fitR=-sum(n1(1:u(j,1)));
for jlevel=2:nd
w2=sum(probR(u(j,jlevel-1)+1:u(j,jlevel)));

n2=(probR/w2).*log(probR/w2);
fitR=fitR-sum(n2(u(j,jlevel-1)+1:u(j,jlevel)));

end
we=sum(probR(u(j,nd)+1:255));
ne=(probR/we).*log(probR/we);
fitR=fitR-sum(ne(u(j,nd)+1:255));

fnew=fitR;
