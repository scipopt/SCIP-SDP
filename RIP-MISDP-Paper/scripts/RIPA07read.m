function [X, rk, ubdelta] = RIPA07read(n, tol, side, file)
%liest SCIP-SDP-Lösung 'file' und liefert Matrix X sowie Rang (mit Tolerenz
%tol) sowie obere Schranke für delta (für
%linke Seite der RIP falls side='l', rechte Seite falls side='r')
%kompatibel mit RIPSDPA und Aspr07SDPA
fid = fopen(file, 'r');
B=zeros(3, n^2+n);
a=fscanf(fid, 'objective value: %f ');
B=fscanf(fid, 'X_%f %f (obj:%f) ');
fclose(fid);
X=zeros(n,n);
for(l=0:1:(length(B)/3)-1)
    [i,j]=MatrixPos(n, B(3*l+1));
    if i<=n && j<=n && i>=0 && j>=0
        X(i,j)=B(3*l+2);
        X(j,i)=B(3*l+2);
    end
end
rk = rank(X, tol)
if side=='l'
    ubdelta = -1*a + 1
end
if side=='r'
    ubdelta = -1*a - 1 %%da Standardform minimiert muss ZF-Wert noch mit -1 multipliziert werden um "echten" Wert zu erhalten
end
end

