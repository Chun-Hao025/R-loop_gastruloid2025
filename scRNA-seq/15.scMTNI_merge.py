import os
import pandas as pd
import numpy as np

###run the merge output script to combine result of split gene chunk after each random sample + beta1/beta2 combination
idx=range(0, 49) ### this depends on how many chunk of gene set you have (check the ogid folder)
ran=pd.read_table('/path/to/GRN/scMTNI/ID.txt', header=None)[0].tolist() ###ID.txt contains the random sample folder name for each set of scMTNI trials.
clusters=["cluster0", "cluster1", "cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9","cluster10",
          "cluster11","cluster12","cluster13","cluster14","cluster15"]
indir='/path/to/GRN/scMTNI/scMTNI_control'
for ran in ran: 
    for cluster in clusters:
        cluster_path=os.path.join(indir, ran)
        merged=[]
        for i in idx:
            file_path=os.path.join(indir, f'{ran}/ogid{i}/{cluster}/fold0/var_mb_pw_k50.txt')
            if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
                print(f"{file_path} does not exist. Skipping.")
                continue  # Skip to the next file if it doesn't exist
            tmp=pd.read_table(file_path, header=None)
            merged.append(tmp)
        if merged:
            dfmerge = pd.concat(merged, axis=0)
            fname = os.path.join(cluster_path, f'{cluster}_{ran}_var_mb_pw_k50_N.txt') ## change the name as the indicator of trial from a specific beta1/beta2 combination
            dfmerge.to_csv(fname, sep='\t', mode='w', index=None, header=False, float_format='%g')
        else:
            print(f"No files found for {ran} {cluster}. Skipping save.")