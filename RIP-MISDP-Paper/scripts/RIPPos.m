function [ p ] = RIPPos(n,i,j)
% liefert Variablennummer für (i,j)-Eintrag für Dimension n, j>=i
p=(n+1)*n/2-(n-i+2)*(n-i+1)/2+j-i+1;  
end

