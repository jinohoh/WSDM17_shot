# Program : do_tuckernn.sh
# Description : Demo HaTen2-TuckerNN

hadoop fs -put test.tensor .
./run_tuckernn.sh 5:4:3 2:2:2 1 1 test.tensor tuckernn_out tuckernn_log 
