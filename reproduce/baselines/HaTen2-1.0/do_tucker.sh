# Program : do_tucker.sh
# Description : Demo HaTen2-Tucker

#hadoop fs -put test.tensor .
#./run_tucker.sh 5:4:3 2:2:2 4 1 test.tensor tucker_test/tucker_out tucker_test/tucker_log 
#hadoop fs -put test.tensor ./test
#./run_tucker.sh 10000:10000:10000:10000 8:8:8:8 4 1 syn4_d1E+04_n1E+05_ot1 tucker_test/tucker_out tucker_test/tucker_log 
./run_tucker.sh 10000:10000:10000:10000 8:8:8:8 4 1 syn4_d1E+04_n1E+05_ot1 tucker_test/tucker_out tucker_test/tucker_log 
