rm testPred/*.db
python -u predictFromFasta.py ~/peptideList.txt fullTrainedMs2Model.h5 testPred/reg.db 0 100000 &> testPred/reg.txt
python -u predictFromFastaThreading.py ~/peptideList.txt fullTrainedMs2Model.h5 testPred/thread.db 0 100000 &>testPred/thread.txt
