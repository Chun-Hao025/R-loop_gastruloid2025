import os
import pandas as pd
import argparse
import gc

def main(args):
    indir = args.indir + '/'
    frac = args.fraction
    filelist = pd.read_table(args.filelist, header=None)
    ##infiles = filelist[1].values
    cells = filelist[0].values
    infiles = [os.path.join(indir, f'{d}_count.tsv') for d in cells] 
    start_seed = args.start_seed

    for seed in range(start_seed, start_seed + args.nseed):
        merged=[]    
        outdir = os.path.join(indir, f'subsample/randseed{seed}/')
        os.makedirs(outdir, exist_ok=True)
        
        for i, cell in enumerate(cells):
            print(f"Processing {infiles[i]} with seed {seed}")
            data = pd.read_table(infiles[i], sep='\t', index_col=0)
            data1 = data.transpose()
            data2 = data1.sample(n=int(data1.shape[0] * frac), replace=False, random_state=seed)
            data3 = data2.transpose()
            
            # Save individual cell data
            outname = os.path.join(outdir, f'{cell}.table')
            data3.to_csv(outname, index=True, index_label='Gene', header=True, sep='\t', float_format='%g')
            merged.append(data3)

            # Free memory by deleting variables that are no longer needed
            del data, data1, data2, data3
            gc.collect()  # Force garbage collection


        print(f"Merged shape for seed {seed}")
        
        # Free memory for the merged DataFrame
        del merged
        gc.collect()  # Force garbage collection

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(

        description=__doc__,

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--filelist', help='filelist', type=str, default='')
    parser.add_argument('--nseed', help='number of subsamples', type=int, default=50)   
    parser.add_argument('--indir', help='indir', type=str, default='')
    parser.add_argument('--fraction', help='fraction of samples', type=float, default=0.5)
    parser.add_argument('--start_seed', help='fraction of samples', type=int, default=1)

    args = parser.parse_args()
    main(args)



