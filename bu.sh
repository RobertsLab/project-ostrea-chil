rsync -av ../project-ostrea-chil \
--exclude='*.sam' \
--exclude='tmp*' \
--exclude='*C_to_T.fastq' \
--exclude='*G_to_A.fastq' \
--exclude='Non_CpG_context*' \
--exclude='.*' --exclude='*/.*' \
sr320@gannet.fish.washington.edu:/volume1/v1_web/owlshell/bu-github/
