function [] = RIPCBFdual(A, k, side, file, Rank, socp, bounds)
% schreibt SDP-File für ganzzahlige RIP-SDP-Relaxierung in dualer Form
% (mit Skalarvariablen) für Matrix A, Ordnung k, schreibt in 'file' 
% side ='l' für linke Seite/alpha_k, side='r' für rechte Seite/beta_k
% Rank = 1 für zusätzliche Rang-NB, sonst Rank = 0
% socp = 1 für zusätzliche gültige SOCP-Ungleichung von Li/Xie, sonst = 0
% bounds = 1 falls für nichtnegative Matrizen A die Schranke 0 <= X_{ij}
% statt -z_j <= X_{ij} für alle i,j benutzt werden soll, sonst = 0
% ACHTUNG: Schreibt untere Dreiecksmatrizen!

    fid = fopen(file, 'w');
    m=length(A(:,1));
    n=length(A(1,:));
    
    % check if A is entrywise nonnegative:
    if bounds == 1 && ~all(A >= 0, 'all')
        bounds = 0;
        fprintf("Setting bounds = 0, since matrix A is not nonnegative!\n");
    elseif bounds == 1 && ~side == 'r'
        bounds = 0;
        fprintf("Setting bounds = 0, since  only works for right side of RIP!\n");
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

    % add scalar variables z and X_{ij}
    fprintf(fid, "VAR\n");
    fprintf(fid, "%d 2\n", n+0.5*n*(n+1));
    fprintf(fid, "L+ %d\n", n);
    if bounds == 0
        fprintf(fid, "F %d\n", 0.5*n*(n+1));
    elseif bounds == 1
        fprintf(fid, "L+ %d\n", 0.5*n*(n+1));
    else
        error("Error: Option <%s> for parameter bounds not valid!\n", bounds);
    end
    fprintf(fid, "\n");


    % objective
    fprintf(fid, "OBJACOORD\n");
    fprintf(fid, "%d\n", 0.5*n*(n+1));
    cnt = n;
    for i = 0:n-1
        for j = 0:i
            if i ~= j
                fprintf(fid, "%d %.15g\n", 0.5*i*(i+1)+j, B(i+1,j+1)+B(j+1,i+1));
                cnt = cnt + 1;
            else
                fprintf(fid, "%d %.15g\n", 0.5*i*(i+1)+j, B(i+1,j+1));
            end
            % check that matrix is symmetric
            if ( abs(B(i+1,j+1) - B(j+1,i+1)) > 1e-6 )
                error("Error: B matrix not symmetric!\n");
            end
        end
    end
    fprintf(fid, "\n");

    % all scalar variables z are binary
    fprintf(fid, "INT\n");
    fprintf(fid, "%d\n", n);
    for j = 0:n-1
        fprintf(fid, "%d\n", j);
    end
    fprintf(fid, "\n");

    % add psd constraints
    fprintf(fid, "PSDCON\n");
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


    % add rank1-constraint
    if ( Rank == 1 )
        fprintf(fid, "PSDCONRANK1\n");
        fprintf(fid, "1\n");
        fprintf(fid, "0\n");
        fprintf(fid, "\n");
    end


    % ------------------------------------------
    % constraints
    fprintf(fid, "CON\n");
    ncons = n*(n+1)+n+2 - bounds*0.5*n*(n+1);
    fprintf(fid, "%d 5\n", ncons);
    fprintf(fid, "L+ %d\n", n);         % -z_j + 1 >= 0
    if bounds == 0
        fprintf(fid, "L- %d\n", 0.5*n*(n+1));  % -z_j - X_{ij} <= 0
    end
    fprintf(fid, "L+ %d\n", 0.5*n*(n+1));       % z_j - X_{ij} >= 0
    fprintf(fid, "L= 1\n");                    % \sum_j X_{jj} == 1
    fprintf(fid, "L- 1\n");                    % \sum_j z_j <= k
    fprintf(fid, "\n");
    
    
    % write ACOORD
    fprintf(fid, "ACOORD\n");
    fprintf(fid, "%d\n", 2*n*(n+1)+3*n - bounds*n*(n+1));
    cnt = 0;       % counts number of specified entries
    conscnt = 0;   % counts number of specified constraints 

    % add upper bounds on z
    for j = 0:n-1
        fprintf(fid, "%d %d -1.0\n", conscnt, j);
        cnt = cnt + 1;
        conscnt = conscnt + 1;
    end

    % add coupling constraints -z_j - X_{ij} <= 0 (only if bounds = 0)
    if bounds == 0
        for i = 0:n-1
            for j = 0:i
                if ( j ~= i )
                    fprintf(fid, "%d %d -0.5\n", conscnt, j);
                else
                    fprintf(fid, "%d %d -1.0\n", conscnt, j);
                end
                fprintf(fid, "%d %d -1.0\n", conscnt, n+0.5*i*(i+1)+j);
                cnt = cnt + 2;
                conscnt = conscnt + 1;
            end   
        end
    end
    
    % add coupling constraints z_j - X_{ij} >= 0
    for i = 0:n-1
        for j = 0:i
            if ( j ~= i )
                fprintf(fid, "%d %d 0.5\n", conscnt, j);
            else
                fprintf(fid, "%d %d 1.0\n", conscnt, j);
            end
            fprintf(fid, "%d %d -1.0\n", conscnt, n+0.5*i*(i+1)+j);
            cnt = cnt + 2;
            conscnt = conscnt + 1;
        end
    end

    % trace constraint
    for i = 0:n-1
        fprintf(fid, "%d %d 1.0\n", conscnt, n+0.5*i*(i+1)+i);
        cnt = cnt +1;
    end
    conscnt = conscnt + 1;

    % sparsity constraint
    for i = 0:n-1   
        fprintf(fid, "%d %d 1.0\n", conscnt, i);
        cnt = cnt + 1;
    end
    conscnt = conscnt + 1;

    fprintf(fid, "\n");
    if (cnt ~= 2*n*(n+1)+3*n - bounds*n*(n+1) || conscnt ~= ncons)
        fprintf("cnt is = %d, should be = %d\n",cnt, 2*n*(n+1)+3*n - bounds*n*(n+1));
        fprintf("conscnt = %d, should be = %d\n",conscnt, ncons);
        error("Error: Something went wrong when writing ACOORD!\n");
    end

    % write BCOORD
    fprintf(fid, "BCOORD\n");
    fprintf(fid, "%d\n", n+2);
    conscnt = 0;
    cnt = 0; 

    % z < = 1
    for j = 0:n-1
        fprintf(fid, "%d 1.0\n", cnt);
        cnt = cnt + 1;
        conscnt = conscnt + 1;
    end

    % nothing to add for coupling constraints -z_j - X_{ij} <= 0 (only if
    % bounds = 0)
    if bounds == 0
        for i = 0:n-1
            for j = 0:i
                conscnt = conscnt + 1;
            end
        end
    end

    
    % nothing to add for coupling constraints z_j - X_{ij} >= 0
    for i = 0:n-1
        for j = 0:i
            conscnt = conscnt + 1;
        end
    end

    % trace constraint
    fprintf(fid, "%d -1.0\n", cnt);
    cnt = cnt + 1;
    conscnt = conscnt + 1;

    % sparsity constraint
    fprintf(fid, "%d %d\n", cnt, -k);
    cnt = cnt + 1;
    conscnt = conscnt + 1;

    fprintf(fid, "\n");

    if ( conscnt ~= ncons || cnt ~= n+2 )
        fprintf("conscnt is = %d, should be %d\n",conscnt, ncons);
        fprintf("cnt = %d, should be %d\n",cnt, n+2);
        error("Error: Something went wrong when writing BCOORD!\n");
    end


    % write HCOORD
    fprintf(fid, "HCOORD\n");
    if socp == 1
        fprintf(fid, "%d\n", 0.5*n*(n+1) + 4*n+n^2);
    else
        fprintf(fid, "%d\n", 0.5*n*(n+1));
    end

    cnt = 0;       % counts number of specified entries
    for i = 0:n-1
        for j = 0:i
            fprintf(fid, "0 %d %d %d 1.0\n", n+0.5*i*(i+1)+j,i,j);
            cnt = cnt + 1;
        end
    end

    if socp == 1
        socpcnt = 0;
        for i = 0:n-1
            fprintf(fid, "%d %d 0 0 0.5\n", i+1,i);
            fprintf(fid, "%d %d %d 0 -0.5\n", i+1,i,n+1);
            fprintf(fid, "%d %d 0 0 0.5\n", i+1,n+0.5*i*(i+1)+i);
            fprintf(fid, "%d %d %d 0 0.5\n", i+1,n+0.5*i*(i+1)+i,n+1);
            socpcnt = socpcnt + 4;
            for j = 0:n-1
                if j <= i
                    fprintf(fid, "%d %d %d 0 1.0\n", i+1,n+0.5*i*(i+1)+j,j+1);
                else
                    fprintf(fid, "%d %d %d 0 1.0\n", i+1,n+0.5*j*(j+1)+i,j+1);
                end
                socpcnt = socpcnt + 1;
            end
        end
    else
        socpcnt = 4*n+n^2;
    end

    if cnt + socpcnt ~= 0.5*n*(n+1) + 4*n+n^2
        error("Error: Something went wrong when writing HCOORD!\n");
    end

    % write DCOORD
    if socp == 1
        fprintf(fid, "DCOORD\n");
        fprintf(fid, "%d\n", n*(n+1));
        cnt = 0;       % counts number of specified entries
        for i = 0:n-1
            for j = 0:n
                fprintf(fid, "%d %d %d 1.0\n", i+1,j+1,j+1);
                cnt = cnt + 1;
            end
        end
        if cnt ~= n*(n+1)
            error("Error: Something went wrong when writing DCOORD!\n");
        end
    end


    % close file
    fclose(fid);    
end
