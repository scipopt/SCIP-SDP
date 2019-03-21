function [ A ] = SDPRelaxCont( m,n,k,option,p,r,MatrixFile, RIPlFile, RIPrFile)
%liefert sdpa-Files für alle SDP-Relaxierungen mit Dateiname =|= ' '
%Parameter p,r nur für DeVore Matrix benötigt, p ebenfalls für Bandmatrix
%option: 'norm': N(0,1)
%'devo': deterministische p² x p^(r+1) DeVore-Matrix die RIP für
%delta=(k-1)r/p erfüllt, m und n irrelevant
%'wish': N(0,1/m)
%'bern': +/- 1/sqrt(m) gleichverteilt
%'0+-1': P(+sqrt(3/m))=1/6, P(0)=2/3, P(-sqrt(3/m))=1/6
%'band': mxm Bandmatrix mit Banbreite p (muss ungerade sein) und gleichverteilten Einträgen in {0,1}
%'rnk1': mxm Matrix aa^T mit Rang 1, a ist N(0,1)-verteilt
%'bina': mxn Binärmatrix, Einträge 0/1-gleichverteilt
fid = fopen(MatrixFile, 'w');
if option=='norm'
    T=[m n k];
    S=['randomization = N(0,1), m= %-4.0g , n= %-4.0g , k= %-4.0g \n'];
    A=zeros(m,n);
    for(i=1:1:m)
        for(j=1:1:n)
            A(i,j)=random('norm',0,1);
            T=[T A(i,j)];
            S= [S '%-4.6g '];
        end
        S=[S '\n'];
    end
end
if option=='devo'
    T=[p r (k-1)*r/p m n k];
    S=['deterministic DeVore Matrix for p=%-4.0g, r=%-4.0g, delta=%-4.3g, m= %-4.0g , n= %-4.0g , k= %-4.0g \n'];
    A=zeros(p^2, p^(r+1));
    f=zeros(1,r+1);
    for l=1:1:p^(r+1)
        for pos=0:1:r
            if mod(l, p^(r-pos+1)) == 0
                f(1,pos+1)=ceil(l/p^(r-pos))-1;
            else
                f(1,pos+1)=ceil(mod(l, p^(r-pos+1))/p^(r-pos))-1; %% berechnet Koeffizienten des Polynoms
            end
        end
        for x=0:1:p-1
            for y=0:1:p-1
                row=1+x*p+y; %% berechnet Position in lexikographischer Ordnung für (x,y)
                sum=0;
                for pos=0:1:r
                    sum=sum+f(1,pos+1)*x^pos; %% Berchnet Funktionswert des Polynoms
                end
                A(row,l)=(y==mod(sum,p)); %% Matrix Eintrag ist 1 falls y=polynomial(x) [in F_p], sonst 0
            end
        end
    end
    A=(1/sqrt(p))*A;
    for i=1:1:length(A(:,1))
        for j=1:1:length(A(1,:))
            T=[T A(i,j)];
            S= [S '%-4.6g '];
        end
        S=[S '\n'];
    end
end
if option == 'wish'
    T=[m n k];
    S=['randomization = N(0,1/m), m= %-4.0g , n= %-4.0g , k= %-4.0g \n'];
    A=zeros(m,n);
    for(i=1:1:m)
        for(j=1:1:n)
            A(i,j)=random('norm',0,sqrt(1/m));
            T=[T A(i,j)];
            S= [S '%-4.6g '];
        end
        S=[S '\n'];
    end
end
if option == 'bern'
    rand=0;
    T=[m n k];
    S=['randomization = uniformly in +/- 1/sqrt(m), m= %-4.0g , n= %-4.0g , k= %-4.0g \n'];
    A=zeros(m,n);
    for(i=1:1:m)
        for(j=1:1:n)
            rand=random('unid',2);
            if rand==1
                A(i,j)=-1*(1/sqrt(m));
            else
                A(i,j)=1/sqrt(m);  
            end
            T=[T A(i,j)];
            S= [S '%-4.6g '];
        end
        S=[S '\n'];
    end
end
if option == '0+-1'
    rand=0;
    T=[m n k];
    S=['randomization = P(+sqrt(3/m))=1/6, P(0)=2/3, P(-sqrt(3/m))=1/6, m= %-4.0g , n= %-4.0g , k= %-4.0g \n'];
    A=zeros(m,n);
    for(i=1:1:m)
        for(j=1:1:n)
            rand=random('unid',6);
            if rand==1
                A(i,j)=sqrt(3/m);
            end
            if rand==6
                A(i,j)=-1*(sqrt(3/m));
            end    
            T=[T A(i,j)];
            S= [S '%-4.6g '];
        end
        S=[S '\n'];
    end
end
if option == 'band'
    rand=0;
    T=[m p k];
    S=['mxm band matrix with entries uniformly in {0,1}, m= %-4.0g , bandwith= %-4.0g , k= %-4.0g \n'];
    A=zeros(m,m);
    for(i=1:1:m)
        for(j=1:1:m)
            rand=random('unid',2);
            if rand==2 && j<=i+floor(p/2) && j >= i-floor(p/2)
                A(i,j)=1;
            end    
            T=[T A(i,j)];
            S= [S '%-4.6g '];
        end
        S=[S '\n'];
    end
end
if option== 'rnk1'
    T=[m k];
    S=['rank 1 Matrix aa^T, randomization for a = N(0,1), m= %-4.0g , k= %-4.0g \n'];
    a=zeros(m,1);
    for i=1:1:m
        a(i)=random('norm',0,1);
    end
    A=a*transpose(a);
    for(i=1:1:m)
        for(j=1:1:n)
            T=[T A(i,j)];
            S= [S '%-4.6g '];
        end
        S=[S '\n'];
    end
end
if option== 'bina'
    rand=0;
    T=[m n k];
    S=['randomization = uniformly in 0/1, m= %-4.0g , n= %-4.0g , k= %-4.0g \n'];
    A=zeros(m,n);
    for(i=1:1:m)
        for(j=1:1:n)
            rand=random('unid',2);
            if rand==1
                A(i,j)=1;
            else
                A(i,j)=0;  
            end
            T=[T A(i,j)];
            S= [S '%-4.6g '];
        end
        S=[S '\n'];
    end
end
fprintf(fid, S, T);
fclose(fid);
if RIPlFile ~= ' '
    RIPSDPAcont(A,k,'l',RIPlFile,0);
end
if RIPrFile ~= ' '
    RIPSDPAcont(A,k,'r',RIPrFile,0);
end
% if Aspr07lFile ~= ' '
%     Aspr07SDPA(A,k,'l',Aspr07lFile,0);
% end
% if Aspr07rFile ~= ' '
%     Aspr07SDPA(A,k,'r',Aspr07rFile,0);
% end
% if nullspaceContFile ~= ' '
%     nullspaceSDPAdualized(A,k,nullspaceContFile,0,0); 
% end
% if RIPlZPLFile ~= ' '
%     RIPAZimpl(A,k,'l',RIPlZPLFile);
% end
% if RIPrZPLFile ~= ' '
%     RIPAZimpl(A,k,'r',RIPrZPLFile);
% end
% if nullspaceZPLFile ~= ' '
%     NSPAZimpl(A,k,nullspaceZPLFile);
% end
end

