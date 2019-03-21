function [ i j ] = MatrixPos(n, l)
%liefert die Position (i,j) für Variable l einer n x n Matrix
l=l+1; %SCIP-Nummerierung startet bei 0 statt 1
i=-1;
j=-1;
for p=1:1:n
    if l <= n-(p-1)
        i=p;
        j=(p-1)+l; %% p-1 Einträge vor erstem aufgeführten Eintrag
        break
    else
        l = l-(n-(p-1));
    end
end

end

