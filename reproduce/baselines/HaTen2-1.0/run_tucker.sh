# Program : run_tucker.sh
# Description : Run HaTen2-Tucker

which hadoop > /dev/null
status=$?
if test $status -ne 0 ; then
	echo ""
	echo "Hadoop is not installed in the system."
	echo "Please install Hadoop and make sure the hadoop binary is accessible."
	exit 127
fi

if [ $# -ne 7 ]; then
	 echo 1>&2 "Usage: $0 [dim_1:..:dim_N (tensor)] [dim_1:...:dim_N (core tensor)] [# of reducers] [max iteration] [tensor path] [output path] [log file]"
	 exit 127
fi

hadoop jar HaTen2.jar kr.ac.snu.datalab.HaTen2.Tucker $1 $2 $3 $4 $5 $6 $7
rm lanczos.* ritz.*
