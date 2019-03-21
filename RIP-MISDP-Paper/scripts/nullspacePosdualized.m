function [ l ] = nullspacePosdualized(p,n,i,j,A)
% liefert Variablennummer für (i,j)-Eintrag von A für Dimension n+p, j>=i
% A <- 'X', 'Y', 'Z', '0'
% negativer Output falls fehlerhafter Input (Details siehe unten)
N=n+p;
if A=='X'
    l=(N+1)*N/2-(N-i+2)*(N-i+1)/2+j-i+1;
end
if A=='Z' %%liefert Position von i,j-Eintrag von Z^T
    jj=i;
    ii=j; %%Variable stammt aus Z^T
     l=(N+1)*N/2-(N-ii+2)*(N-ii+1)/2+jj-ii+1+p; 
end
if A== 'Y'
     l=(N+1)*N/2-(N-(i+p)+2)*(N-(i+p)+1)/2+j-i+1;
end
if ((i>p | j>p) && A=='X') 
    l=-1; %%%irgendwas schiefgelaufen
end
if ((i>n | j>n) && A=='Y')
    l=-2; %%%irgendwas schiefgelaufen
end
if ((i>n | j>p) && A=='Z')
    l=-3; %%%irgendwas schiefgelaufen
end
if (i>j && A ~= 'Z')
    l=-4; %%%irgendwas schiefgelaufen
end
if A== '0' %%%gesamte Matrix
    l=(N+1)*N/2-(N-i+2)*(N-i+1)/2+j-i+1;  
end
end

