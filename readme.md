I.	How to use the script.

1) Place “start_csv.csv” (see attached example file “start_csv.csv”) and “dirty_seqs.r” files into the same directory.
2) Open the script in “R” (it was made under R v. 4.3.2).
3) Start the script.

II.	Some points.

1) Be sure that all necessary libraries are installed and working properly.
2) You should change “gap.end” and “gap.mid” parameters according to your inner reasons.
3) Remove the “#” if you want to test all models and to use the best fitted model. Do not forget to replace the best fitted model downstream.

III.	If you want to start the script under Windows OS.

1) Place “start_csv.csv”, “dirty_seqs.r”, and “start_dirty_seqs.bat” files into the same directory.
2) Change “start_dirty_seqs.bat” file according to your local OS (you must replace the path to directory containing “Rscript.exe”). 
NB! You'll probably have to correct the script according to the command prompt of your local operating system. 
If it done correctly, working directory will be set automatically.
3) Start the “start_dirty_seqs.bat” with double click.

IV.	The script will generate four files. It could be useful during the processing of your data:
  - “final_output.fasta” – *.fasta file after downloading the sequences from GenBank NCBI.
  - “myMuscleAlignment.fasta” – aligned sequences.
  - “trimmedMuscleAlignment.fasta” – trimmed sequences.
  - “ML_tree.tre” – Maximum Likelihood tree with rapid bootstraps based on trimmed alignment.

Feel free to contact me: svshchenkov@yandex.ru (Sergei).
