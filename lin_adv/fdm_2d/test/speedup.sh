# Script to perform speedup test
if [ $# -lt 1 ]
then
   echo "Usage: Specify list of processors"
   echo "       sh ./speedup.sh \"1 2 3 4 5 6 7 8 9 10\""
   exit 1
fi

procs=$1

rm -f *.plt *.h5 times.txt
for p in $procs
do
   echo "Running with $p processors"
   mpirun -np $p ../src/exe -log_view > log_$p.txt || echo "mpirun failed"
   t=`grep "Iter runtime min" log_$p.txt | awk '{print $5}'`
   echo "$p  $t" >> times.txt
   rm -rf *.plt *.h5
done

