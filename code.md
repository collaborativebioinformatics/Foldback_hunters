## Simulating HG002 reads

Code 
```
srun \
    --job-name=pbsim_chr6 \
    --time=03:00:00 \
    --cpus-per-task=4 \
    --mem=16G \
    --output=pbsim_chr6.log \
    pbsim \
        --strategy wgs \
        --method errhmm \
        --errhmm ERRHMM-ONT-HQ.model \
        --depth 50 \
        --length-mean 15000 \
        --length-sd 8000 \
        --length-min 1000 \
        --length-max 200000 \
        --genome HG002_chr6_both.fasta \
        --prefix chr6_sim
```
