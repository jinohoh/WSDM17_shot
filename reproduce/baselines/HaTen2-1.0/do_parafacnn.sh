# Program : do_parafacnn.sh
# Description : Demo HaTen2-PARAFACNN

hadoop fs -put test.tensor .
./run_parafacnn.sh 5:4:3 2 1 1 test.tensor parafacnn_out parafacnn_log 
