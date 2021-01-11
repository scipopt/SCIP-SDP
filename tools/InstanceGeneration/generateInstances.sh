# generate random instances according to the paper from Kobayashi and
# Takano (2020): A branch-and-cut algorithm for solving mixed-integer
# semidefinite optimization problems

for n in {15,30,60};
do
    perl generateRandomMISDP.prl $n 30 30 1.0 5;
    # mv randomMISDP_${n}_30_30_1.0.cbf randomMISDP_${n}_30_30_1.0_${cnt}.cbf;

    perl generateRandomMISDP.prl 30 $n 30 1.0 5;
    # mv randomMISDP_30_${n}_30_1.0.cbf randomMISDP_30_${n}_30_1.0_${cnt}.cbf;

    perl generateRandomMISDP.prl 30 30 $n 1.0 5;
    # mv randomMISDP_30_30_${n}_1.0.cbf randomMISDP_30_30_${n}_1.0_${cnt}.cbf;
done

for a in {0.1,10.0};
do
    perl generateRandomMISDP.prl 30 30 30 $a 5;
    # mv randomMISDP_30_30_30_${a}.cbf randomMISDP_30_30_30_${a}_${cnt}.cbf;
done

for n in {45,60};
do
    perl generateRandomMISDP.prl $n $n $n 1.0 5;
    # mv randomMISDP_${n}_${n}_${n}_1.0.cbf randomMISDP_${n}_${n}_${n}_1.0_${cnt}.cbf;
done

    
for x in *1.0*;
do
    mv $x ${x/1.0/1}
done

for x in *10.0*;
do
    mv $x ${x/10.0/10}
done

for x in *0.1*;
do
    mv $x ${x/0.1/0_1}
done

for x in *.cbf;
do
    gzip $x;
done
