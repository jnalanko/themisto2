set -xue

mkdir -p export
mkdir -p fof
mkdir -p index
mkdir -p logs
mkdir -p sbwt  
mkdir -p seqs  
mkdir -p temp
mkdir -p unitigs

ls seqs/ | xargs --verbose -P 16 -I {} ggcat build -p -s 1 -m 2 --temp-dir temp -j 2 -k 31 seqs/{} -o unitigs/{}.unitigs.fna

SEQ_COUNT=$(find seqs -maxdepth 1 -type f | wc -l)
HALF_COUNT=$(( SEQ_COUNT / 2 ))

find seqs -type f | grep "\.fna" | sort > fof/seqs.txt # Might be gzipped

find unitigs/ -type f | grep ".fna$" | sort > fof/unitigs.txt # Are not gzipped because they come from ggcat

find unitigs/ -type f | grep ".fna$" | sort | head -n $HALF_COUNT > fof/unitigs-half1.txt
find unitigs/ -type f | grep ".fna$" | sort | tail -n $HALF_COUNT > fof/unitigs-half2.txt

echo "index/half1.thm2" > fof/merge.txt
echo "index/half2.thm2" >> fof/merge.txt
