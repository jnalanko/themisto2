set -xueo pipefail

/usr/bin/time -v sbwt build --input-list fof/unitigs.txt -v -o sbwt/unitigs --temp-dir temp -l -k 31 -m 50 -t 32 -r --in-memory 2>&1 | tee logs/sbwt.log
/usr/bin/time -v sbwt build --input-list fof/unitigs-half1.txt -v -o sbwt/unitigs-half1 --temp-dir temp -l -k 31 -m 50 -t 32 -r --in-memory 2>&1 | tee logs/sbwt-half1.log
/usr/bin/time -v sbwt build --input-list fof/unitigs-half2.txt -v -o sbwt/unitigs-half2 --temp-dir temp -l -k 31 -m 50 -t 32 -r --in-memory 2>&1 | tee logs/sbwt-half2.log

/usr/bin/time -v themisto2 build -s sbwt/unitigs.sbwt -i fof/seqs.txt -o index/seqs.thm2 --temp-dir temp -k 31 -d 30 -t 32 --index-type sparse-dense 2>&1 | tee logs/seqs.log

/usr/bin/time -v themisto2 build -s sbwt/unitigs.sbwt -i fof/unitigs.txt -o index/unitigs.thm2 --temp-dir temp -k 31 -d 30 -t 32 --index-type sparse-dense 2>&1 | tee logs/unitigs.log

/usr/bin/time -v themisto2 build -s sbwt/unitigs-half1.sbwt -i fof/unitigs-half1.txt -o index/half1.thm2 --temp-dir temp -k 31 -d 30 -t 32 --index-type sparse-dense 2>&1 | tee logs/half1.log

/usr/bin/time -v themisto2 build -s sbwt/unitigs-half2.sbwt -i fof/unitigs-half2.txt -o index/half2.thm2 --temp-dir temp -k 31 -d 30 -t 32 --index-type sparse-dense 2>&1 | tee logs/half2.log

/usr/bin/time -v themisto2 merge --index-file-list fof/merge.txt -o index/merge.thm2 --temp-dir temp -d 30 -t 32 2>&1 | tee logs/merge.log

# Export
/usr/bin/time -v themisto2 export -i index/seqs.thm2 -o export/seqs -t 32 2>&1 | tee logs/half2.log
/usr/bin/time -v themisto2 export -i index/unitigs.thm2 -o export/seqs -t 32 2>&1 | tee logs/half2.log
/usr/bin/time -v themisto2 export -i index/merge.thm2 -o export/seqs -t 32 2>&1 | tee logs/half2.log

