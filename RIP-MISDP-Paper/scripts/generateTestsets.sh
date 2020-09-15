if [ ! -f MISDP_backup/RIPSDPA.solu ]; then
    echo "ERROR: File /MISDP_backup/RIPSDPA.solu is missing!";
    return
fi


cd MISDP/

fileSDPA=RIPSDPA.test
for x in `ls *.dat-s`;
do
    printf "MISDP/RIPtest/%s\n" $x >> $fileSDPA;
done;

cnt_type=0;
cnt_instance=0;
testfileAll=RIPCBFall.test;
solufileAll=RIPCBFall.solu;

# options:
arr_rank1=(0 1);
arr_trineq=(0 1);
arr_boundver=(0 1 2);
arr_pd=("p" "d");
arr_socp=(0 1);
arr_strgbnds=(0 1);

# rank1-constraint
for r in "${arr_rank1[@]}";
do
    # trace-inequality
    for t in "${arr_trineq[@]}";
    do
	# boundver
	for b in "${arr_boundver[@]}";
	do
	    # primal/dual
	    for pd in "${arr_pd[@]}";
	    do
		# socp-inequality
		for so in "${arr_socp[@]}";
		do
		    # strong bounds
		    for st in "${arr_strgbnds[@]}";
		    do
			cmd="ls *\$pd\$r\$so\$st\$t\$b.cbf"
			testfile=RIPCBF$pd$r$so$st$t$b.test
			solufile=RIPCBF$pd$r$so$st$t$b.solu
			type=$pd$r$so$st$t$b;
			for x in `eval $cmd`
			do
			    newname="${x/$type.cbf/}"
			    grep $newname ../MISDP_backup/RIPSDPA.solu >> $solufile;
			    grep $newname ../MISDP_backup/RIPSDPA.solu >> $solufileAll;
			    printf "MISDP/RIPtest/%s\n" $x >> $testfile;
			    printf "MISDP/RIPtest/%s\n" $x >> $testfileAll;
			    ((cnt_instance+=1));
			done
			((cnt_type+=1))
		    done
		done
	    done
	done
    done
done

printf "Created overall testset with %d instances and %d testsets for the different RIP formulations.\n" $cnt_instance $cnt_type;

cd ..

