# Program : do_parafac.sh
# Description : Demo HaTen2-PARAFAC

hadoop fs -put test.tensor .
./run_parafac.sh 5:4:3 2 1 1 test.tensor parafac_out parafac_log 
