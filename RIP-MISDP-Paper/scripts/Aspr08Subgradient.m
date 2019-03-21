function [ alphabetak ] = Aspr08Subgradient( A,k,side,maxIt,outputfile )
%berechnet Schranke für alpha_k/beta_k via Subgradienten-Verfahren
%Lagrange-Relaxierung wie in d'Aspremont 2008
%A,k wie dort, side <- {'l','r'}: linke oder rechte Seite der RIP
tic
s=1;
h=1;
rho=0;
n=length(A(1,:));
X=zeros(n,n); 
alpha=max(max(A))+1; % größer als größter Matrixeintrag, +1 zur Sicherheit vs. numerische Instabilität
noimprovements=0;
T=zeros(1,maxIt+1);
totalTime=0;
L=0;
while norm(h)>0.0001 && s <=maxIt && noimprovements <= 9 && totalTime < 3600 % Stopp falls h=0, max It oder 10x nicht verbessert 
    if s==1
        fid=fopen(outputfile, 'w');
        S='next Iteration= %-4.6g\n\n';
        P=[s];
        fprintf(fid, S, P);
        fclose(fid);
    else
        fid=fopen(outputfile, 'a');
        S='\ntime= %-4.6g, rho= %-4.6g Value of Lagrange Relaxation= %-4.6g, next Iteration= %-4.6g, \n\n';
        P=[totalTime rho L s];
        fprintf(fid, S, P);
        fclose(fid);
    end
    stepsize=0.05*(1/sqrt(s));
    T(maxIt+1)=T(maxIt+1)+toc;
    tic;
    [Y,L]=Aspr08LagrangeRelaxSDP(A,k,side,rho,alpha,outputfile);
    T(s)=toc;
    tic
    z=0;
    for i=1:n
       for j=1:n
           z=z+1;
           X(i,j)=Y(z);
       end
    end           
    L
    h=k-sqrt(MatrixCard(X,0.001)); % falls Rank=1 gilt CardX=cardx², sonst bestmögliche Approximation  
    rho= rho - stepsize * h 
    if rho <0
        rho = 0;
    end
    if s==1
        alphabetak=L;
    end
    if s>1 && side=='r'
        alphabetak=min([real(alphabetak), real(L)]); % sedumi gibt oftmals komplexe Zahlen mit sehr kleinem (< 10^-8) Imaginärteil aus, Matlab berechnet dann aber Maximum der Beträge statt Maximum, schlecht bei negativen Werten
    end
    if s>1 && side=='l'
        alphabetak=max([real(alphabetak), real(L)]); 
    end
    if s==1
        alphabetakOld = alphabetak;
    else
        if alphabetakOld == alphabetak
            noimprovements = noimprovements + 1;
        else
            alphabetakOld = alphabetak;
            noimprovements = 0;
        end
    end
    s=s+1
    totalTime=sum(T)
    averageTime=(sum(T)-T(maxIt+1))/(s-1)% s-1 da s am Ende nochmal hochgezählt wird

end
finished=1
T(maxIt+1)=T(maxIt+1)+toc;
totalTime=sum(T)
averageTime=(sum(T)-T(maxIt+1))/(s-1) % s-1 da s am Ende nochmal hochgezählt wird
fid=fopen(outputfile, 'a');
S='\ntotal Iterations %-4.6g, total time= %-4.6g, average time per SDP = %-4.6g, best found bound=%-4.6g\n\n';
P=[s-1 totalTime averageTime alphabetak];
fprintf(fid, S, P);
fclose(fid);
end

