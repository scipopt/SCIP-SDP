function [] = RIPCBFprimal(A, k, side, file, Rank, socp, bounds)
% schreibt SDP-File für ganzzahlige RIP-SDP-Relaxierung in primaler Form
% (mit PSD-Variablen) für Matrix A, Ordnung k, schreibt in 'file' 
% side ='l' für linke Seite/alpha_k, side='r' für rechte Seite/beta_k
% Rank = 1 für zusätzliche Rang-NB, sonst Rank = 0
% socp = 1 für zusätzliche gültige SOCP-Ungleichung von Li/Xie, sonst = 0
% bounds = 1 falls für nichtnegative Matrizen A die Schranke 0 <= X_{ij}
% statt -z_j <= X_{ij} für alle i,j benutzt werden soll, sonst = 0
% ACHTUNG: schreibt untere Dreiecksmatrizen!

    fid = fopen(file, 'w');
    m=length(A(:,1));
    n=length(A(1,:));

    % check if A is entrywise nonnegative:
    if bounds == 1 && ~all(A >= 0, 'all')
        bounds = 0;
        fprintf("Setting bounds = 0, since matrix A is not nonnegative!\n");
    elseif bounds == 1 && ~side == 'r'
        bounds = 0;
        fprintf("Setting bounds = 0, since this only works for right side of RIP!\n");
    end
    
    % SOCP-inequality is only valid for right side of the RIP
    if socp == 1 && ~side == 'r'
        socp = 0;
        fprintf("Setting socp = 0, since this only works for right side of RIP!\n");
    end

    % compute B = A^T A
    B = transpose(A)*A;

    % output file
    fprintf("Writing problem to file <%s> ...\n", file);

    fprintf(fid, "VER\n1\n\n");

    if side == 'l'
        fprintf(fid, "OBJSENSE\nMIN\n\n");
    elseif side == 'r'
        fprintf(fid, "OBJSENSE\nMAX\n\n");
    else
        error("Error: Option <%s> for parameter side not valid!\n", side);
    end

    %% add scalar variables z
    fprintf(fid, "VAR\n");
    fprintf(fid, "%d 1\n", n);
    fprintf(fid, "L+ %d\n", n);
    fprintf(fid, "\n");

    %% add Matrix variables X
    fprintf(fid, "PSDVAR\n");
    if socp == 0
        fprintf(fid, "1\n");
        fprintf(fid, "%d\n", n);
    elseif socp == 1
        fprintf(fid, "%d\n", n+1);
        fprintf(fid, "%d\n", n);
        for j = 0:n-1
            fprintf(fid, "%d\n", n+2);
        end
    else
        error("Error: Option <%s> for parameter socp not valid!\n", socp);
    end
    fprintf(fid, "\n");

    %% objective
    fprintf(fid, "OBJFCOORD\n");
    fprintf(fid, "%d\n", 0.5*n*(n+1));
    for i = 0:n-1
        for j = 0:i
            fprintf(fid, "0 %d %d %.15g\n", i, j, B(i+1,j+1));
            % check that matrix is symmetric
            if ( abs(B(i+1,j+1) - B(j+1,i+1)) > 1e-6 )
                error("Error: B matrix not symmetric!\n");
            end
        end
    end

    fprintf(fid, "\n");

    %% all scalar variables z are binary
    fprintf(fid, "INT\n");
    fprintf(fid, "%d\n", n);
    for j = 0:n-1
        fprintf(fid, "%d\n", j);
    end
    fprintf(fid, "\n");

    %% add rank1-constraint
    if ( Rank == 1 )
        fprintf(fid, "PSDVARRANK1\n");
        fprintf(fid, "1\n");
        fprintf(fid, "0\n");
        fprintf(fid, "\n");
    end

    %% constraints
    fprintf(fid, "CON\n");
    nsocpcons = 0.5*n*(n+3)*(n+2);
    ncons = n*(n+1)+n+2 + socp*nsocpcons;
    ncones = 5 + socp;
    fprintf(fid, "%d %d\n", ncons, ncones);
    if socp == 1
        % SOCP constraint is the first constraint
        fprintf(fid, "L= %d\n", nsocpcons);  % \sum_j X_{ij}^2 <= X_{ii}z_i
    end
    fprintf(fid, "L+ %d\n", n);              % -z_j + 1 >= 0
    fprintf(fid, "L- %d\n", 0.5*n*(n+1));    % -z_j - X_{ij} <= 0
    fprintf(fid, "L+ %d\n", 0.5*n*(n+1));    % z_j - X_{ij} >= 0
    fprintf(fid, "L= 1\n");                  % \sum_j X_{jj} == 1
    fprintf(fid, "L- 1\n");                  % \sum_j z_j <= k
    fprintf(fid, "\n");

    %% write FCOORD
    fprintf(fid, "FCOORD\n");
    nFCOORD = n*(n+1)+n + socp*(n*(3*n+5+0.5*n*(n+1)));
    fprintf(fid, "%d\n", nFCOORD);
    cnt = 0;
    conscnt = 0;
    
    % SOCP constraint
    if socp == 1
        for i = 0:n-1
            for j = 0:n
                fprintf(fid, "%d %d %d %d 1.0\n",conscnt,i+1,j+1,j+1);
                cnt = cnt + 1;
                conscnt = conscnt + 1;
            end
            fprintf(fid, "%d %d 0 0 1.0\n",conscnt,i+1);
            fprintf(fid, "%d 0 %d %d -0.5\n",conscnt,i,i);
            conscnt = conscnt + 1;
            cnt = cnt + 2;
            for j = 0:n-1
                if j == i
                    fprintf(fid, "%d %d %d 0 0.5\n",conscnt,i+1,j+1);
                else
                    fprintf(fid, "%d %d %d 0 1.0\n",conscnt,i+1,j+1); 
                end
                if j <= i
                    fprintf(fid, "%d 0 %d %d -1.0\n",conscnt,i,j);
                else
                    fprintf(fid, "%d 0 %d %d -1.0\n",conscnt,j,i);
                end
                conscnt = conscnt + 1;
                cnt = cnt + 2;
            end
            fprintf(fid, "%d %d %d 0 0.5\n",conscnt,i+1,n+1);
            fprintf(fid, "%d 0 %d %d -0.5\n",conscnt,i,i);
            conscnt = conscnt + 1;
            cnt = cnt + 2;
            for k = 0:n
                for j = 0:k-1
                    fprintf(fid, "%d %d %d %d 1.0\n",conscnt,i+1,k+1,j+1);
                    conscnt = conscnt + 1;
                    cnt = cnt + 1;
                end
            end
        end
    end
    if ( conscnt ~= socp*nsocpcons || cnt ~= socp*(n*(3*n+5+0.5*n*(n+1))) )
        fprintf("conscnt is = %d, should be = %d\n",conscnt, socp*nsocpcons);
        fprintf("cnt is = %d, should be = %d\n",cnt,socp*(n*(3*n+5+0.5*n*(n+1))));
        error("Error while writing FCOORD for SOCP constraint!\n");
    end

    % no FCOORD for -z_j + 1 >= 0
    conscnt = conscnt + n;

    % add coupling constraints -z_j - X_{ij} <= 0
    for i = 0:n-1
        for j = 0:i
            if i ~= j
                fprintf(fid, "%d 0 %d %d -0.5\n", conscnt, i, j);
            else
                fprintf(fid, "%d 0 %d %d -1.0\n", conscnt, i, j);
            end
            cnt = cnt + 1;
            conscnt = conscnt + 1;
        end
    end

    % add coupling constraints z_j - X_{ij} >= 0
    for i = 0:n-1
        for j = 0:i
            if i ~= j
                fprintf(fid, "%d 0 %d %d -0.5\n", conscnt, i, j);
            else
                fprintf(fid, "%d 0 %d %d -1.0\n", conscnt, i, j);
            end
            cnt = cnt + 1;
            conscnt = conscnt + 1;
        end
    end

    % trace constraint
    for i = 0:n-1
        fprintf(fid, "%d 0 %d %d 1.0\n", conscnt, i, i);
        cnt = cnt + 1;
    end
    conscnt = conscnt + 1;

    % no FCOORD for sparsity constraint
    conscnt = conscnt + 1;

    fprintf(fid, "\n");

    if ( conscnt ~= ncons || cnt ~= nFCOORD )
        error("Error: Something went wrong when writing FCOORD!\n");
    end

    %% ACOORD
    fprintf(fid, "ACOORD\n");
    if bounds ~= 0 && bounds ~= 1
        error("Error: Option <%s> for parameter bounds not valid!\n", bounds);
    end
    nACOORD = n*(n+1)+2*n + socp*(2*n) - bounds*(0.5*n*(n+1));
    fprintf(fid, "%d\n", nACOORD);
    cnt = 0;
    conscnt = 0;

    % SOCP constraint
    if socp == 1
        for i = 0:n-1
            pos = i*(0.5*(n+3)*(n+2));
            fprintf(fid, "%d %d -0.5\n",pos+n+1,i);
            cnt = cnt + 1;
            fprintf(fid, "%d %d 0.5\n",pos+2*n+2,i);
            cnt = cnt + 1;
        end        
    end

    conscnt = socp*nsocpcons;
    if cnt ~= socp*(2*n)
        error("Error while writing ACOORD for SOCP constraint!\n");
    end  
    
    % add upper bounds on z
    for j = 0:n-1
        fprintf(fid, "%d %d -1.0\n", conscnt, j);
        conscnt = conscnt + 1;
        cnt = cnt + 1;
    end

    % add coupling constraints -z_j - X_{ij} <= 0 (if bounds = 1, this
    % changes to -X_{ij} <= 0, so that nothing needs to be specified here
    if bounds == 0
        for i = 0:n-1
            for j = 0:i
                if ( i ~= j )
                    fprintf(fid, "%d %d -0.5\n", conscnt, j);
                else
                    fprintf(fid, "%d %d -1.0\n", conscnt, j);
                end
                conscnt = conscnt + 1;
                cnt = cnt + 1;
            end
        end
    else
        conscnt = conscnt + 0.5*n*(n+1);
    end
    
    % add coupling constraints z_j - X_{ij} >= 0
    for i = 0:n-1
        for j = 0:i
            if ( i ~= j )
                fprintf(fid, "%d %d 0.5\n", conscnt, j);
            else
                fprintf(fid, "%d %d 1.0\n", conscnt, j);
            end
            conscnt = conscnt + 1;
            cnt = cnt + 1;
        end
    end

    % no ACOORD for trace constraint
    conscnt = conscnt + 1;

    % sparsity constraint
    for j = 0:n-1   
        fprintf(fid, "%d %d 1.0\n", conscnt, j);
        cnt = cnt + 1;
    end
    conscnt = conscnt + 1;

    fprintf(fid, "\n");
    if ( conscnt ~= ncons || cnt ~= nACOORD )
        error("Error: Something went wrong when writing ACOORD!\n");
    end

    %% write BCOORD
    fprintf(fid, "BCOORD\n");
    nBCOORD = n+2 + socp*(n*(n+1));
    fprintf(fid, "%d\n", nBCOORD);
    conscnt = 0;
    cnt = 0;

    % SOCP constraint
    if socp == 1
        for i = 0:n-1
            for j = 0:n
                pos = i*(0.5*(n+3)*(n+2));
                fprintf(fid, "%d -1.0\n",pos+j);
                cnt = cnt + 1;
            end
        end
    end
    if cnt ~= socp*n*(n+1)
        error("Error while writing BCOORD for SOCP constraint!\n");
    end        
    conscnt = socp*nsocpcons;
    
    % z < = 1
    for j = 0:n-1
        fprintf(fid, "%d 1.0\n", conscnt);
        cnt = cnt + 1;
        conscnt = conscnt + 1;
    end

    % nothing to add for coupling constraints -z_j - X_{ij} <= 0
    for i = 0:n-1
        for j = 0:i
            conscnt = conscnt + 1;
        end
    end

        
    % nothing to add for coupling constraints z_j - X_{ij} >= 0
    for i = 0:n-1
        for j = 0:i
            conscnt = conscnt + 1;
        end
    end


    % trace constraint
    fprintf(fid, "%d -1.0\n", conscnt);
    conscnt = conscnt + 1;
    cnt = cnt + 1;

    % sparsity constraint
    fprintf(fid, "%d %d\n", conscnt, -k);
    conscnt = conscnt + 1;
    cnt = cnt + 1;

    fprintf(fid, "\n");
    if ( conscnt ~= ncons || cnt ~= nBCOORD )
        error("Error: Something went wrong when writing BCOORD!\n");
    end

    %% close file
    fclose(fid);
end
