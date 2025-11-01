%___________________________________________________________________%
% SSNMRA source code demo version 1.0      %                            %
%___________________________________________________________________%
function [NMRbest,fmax,bb]=SSNMRA(n,maxiter,Lb,Ub,fitnessR,probR,d)
disp('SSNMRA is optimizing your problem')
bp=0.05;
breeders=n/5;        % breeder population
counter=0;
% workers=n-breeders;  % worker population
iter=1;
for i=1:n,
    %     Lb+(Ub-Lb).*rand(1,d)
    NMRsolution(i,:)=Lb+(Ub-Lb).*rand(1,d);
end

NMRsolution=sort(NMRsolution);
    NMRsolution=fix(NMRsolution);
        NMRfitness=zeros(n,1);
        for i=1:n
        NMRfitness(i,:)=fobj(NMRsolution(i,:),d,probR);
        end

[fmax,I]=max(NMRfitness);
NMRbest=NMRsolution(I,:);
S=NMRsolution;

Fnew=fmax;
while (iter<=maxiter) % Loop over worker and breeders 
    %For workers phae
  lambda=rand;
    for i=((n/5)+1):n,
        ab=randperm(n);
        c1 = 2*exp(-(4*iter/maxiter)^2); 
        c2=rand();
        c3=rand();
        L=Levy(d);
        r1=rand(); 
        r2=rand(); 
        Fc=2-iter*((2)/maxiter);     
        A1=2*Fc*r1-Fc; 
        C1=2*r2; 
        b=1;             
        ll=(Fc-1)*rand()+1;
       if iter<maxiter/3
            S(i,:)=(NMRsolution(i,:)+L.*(NMRsolution(ab(1),:)-NMRsolution(ab(2),:)));
       elseif iter>maxiter/3 && iter<3*maxiter/4
             if c3<0.5 
                    S(i,:)=NMRbest+c1*((Ub-Lb)*c2+Lb);
             else
                    S(i,:)=NMRbest-c1*((Ub-Lb)*c2+Lb);
             end
       else
            D_alphs=Fc*S(i,:)+A1*((NMRbest-S(i,:)));                   
            X1=D_alphs*exp(b.*ll).*cos(ll.*2*pi)+NMRbest;
            S(i,:)=X1;
       end
       Flag4Ub=S(i,:)>Ub;
        Flag4Lb=S(i,:)<Lb;
        S(i,:)=(S(i,:).*(~(Flag4Ub+Flag4Lb)))+Ub.*Flag4Ub+Lb.*Flag4Lb;
         S(i,:)=sort( S(i,:));
         S(i,:)=fix( S(i,:));
         Fnew=fobj(S(i,:),d,probR);
        
        if (Fnew>=NMRfitness(i)),
            NMRsolution(i,:)=S(i,:);
            NMRfitness(i)=Fnew;
        end
    end
    %For Breeders
    for z=1:breeders;
        if rand>bp
            NMRneighbours=randperm(breeders);
            S(z,:)=(1-lambda).*(S(z,:))+(lambda*(NMRbest-NMRsolution(NMRneighbours(1),:)));

            Flag4Ub=S(z,:)>Ub;
        Flag4Lb=S(z,:)<Lb;
        S(z,:)=(S(z,:).*(~(Flag4Ub+Flag4Lb)))+Ub.*Flag4Ub+Lb.*Flag4Lb;
            S(z,:)=sort( S(z,:));
            S(z,:)=fix( S(z,:));
            
            % Evaluate new solutions
            Fnew=fobj(S(z,:),d,probR);

            if (Fnew>=NMRfitness(z)),
                NMRsolution(z,:)=S(z,:);
                NMRfitness(z)=Fnew;
            end
        end
    end
    if  Fnew == NMRfitness,
        counter=counter+1;
    end

    [fmax,I]=max(NMRfitness);
    NMRbest=NMRsolution(I,:);
    S=NMRsolution;
    bb(iter)=fmax;
    iter=iter+1;
    if counter <=10;
        beta=1;
        sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
        for j=1:n,
            ab=randperm(n);
            a=2-j*((2)/maxiter); % a decreases linearly fron 2 to 0
            s=S(j,:);
            u=randn(size(s))*sigma;
            v=randn(size(s));
            step=u./abs(v).^(1/beta);
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            A1=2*a*r1-a; 
            C1=2*r2;
            D_alpha=abs(C1*NMRbest-S(j,:)); 
            X1=NMRbest-A1*D_alpha;
            r1=rand();
            r2=rand();
            A2=2*a*r1-a; 
            C2=2*r2; 
            D_beta=abs(C2*NMRbest-S(j,:)); 
            X2=NMRbest-A2*D_beta; 
            r1=rand();
            r2=rand();
            A3=2*a*r1-a; 
            C3=2*r2; 
            D_delta=abs(C3*NMRbest-S(j,:)); 
            X3=NMRbest-A3*D_delta; 
            S(j,:)=(X1+X2+X3)/3;

            Flag4Ub=S(j,:)>Ub;
            Flag4Lb=S(j,:)<Lb;
            S(j,:)=(S(j,:).*(~(Flag4Ub+Flag4Lb)))+Ub.*Flag4Ub+Lb.*Flag4Lb;

            S(j,:)=sort( S(j,:));
            S(j,:)=fix( S(j,:));
            % Evaluate new solutions
            Fnew=fobj(S(j,:),d,probR);
                   
             if (Fnew>=NMRfitness(j)),
                NMRsolution(j,:)=S(j,:);
                NMRfitness(j)=Fnew;
            end
        end
    end
end

function L=Levy(d)
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);
L=0.01*step;
