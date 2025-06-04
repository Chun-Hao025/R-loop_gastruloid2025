import os
import pandas as pd
import sys
import numpy as np
import argparse
import re

def genOGid(frames, cells, indir, regdf, splitgene):
    gene = frames[0]
    for i in range(1, len(cells)):
        gene = gene.merge(frames[i], how='outer', on='Gene')
    gene.fillna('None', inplace=True)
    n = gene.shape[0]
    Gene_OGID = 'OG' + np.arange(1, n+1).astype(str).astype('object') + '_1'
    NAME = gene.iloc[:, 1]
    for i in range(2, len(cells) + 1):
        NAME = NAME + ',' + gene.iloc[:, i]
    ogid = pd.DataFrame({'Gene_OGID': Gene_OGID, 'NAME': NAME})
    outname = os.path.join(indir, 'testdata_ogids.txt')
    print('OGID file:', outname)
    ogid.to_csv(outname, index=None, header=True, sep='\t', mode='w', float_format='%g')
    
    if regdf:
        regdata = regdf[0]
        for i in range(1, len(cells)):
            regdata = regdata.merge(regdf[i], how='outer', on='Gene')
        regulators = regdata['Gene']
        tfs = gene[gene['Gene'].isin(regulators)].index + 1
        tfs = pd.DataFrame({'TF': tfs})
        tffilename = os.path.join(indir, 'TFs_OGs.txt')
        print('Regulator OG ID file:', tffilename)
        tfs.to_csv(tffilename, index=None, header=None, sep='\t', mode='w', float_format='%g')

def main(args):
    indir = os.path.join(args.indir, '')
    net_dir = os.path.join(args.net_dir, '')
    outdir = os.path.join(args.outdir, '')
    regfile = args.regfile if args.regfile else None
    splitgene = args.splitgene
    cellsuffix = args.addcellsuffix
    
    files = pd.read_table(args.filelist, header=None)
    cells = files[0].values  
    filelist = [os.path.join(indir, f'{d}.table') for d in cells]  
    
    files[0].to_csv(os.path.join(indir, 'celltype_order.txt'), sep='\n', mode='w', index=None, header=False, float_format='%g')
    
    frames = []
    regdf = []
    
    reg = pd.read_table(regfile, header=None) if regfile else None

    ogid_dir = os.path.join(indir, 'ogids')
    sub = [f for f in os.listdir(f'{indir}ogids/') if os.path.isfile(os.path.join(f'{indir}ogids/', f))]
    number1 = [int(re.search(r'\d+', filename).group()) for filename in sub if re.search(r'\d+', filename)]

    # Loop over gene splits and create individual config files
    for idx, num in zip(sub, number1):
        ogid_outdir = os.path.join(outdir, f'ogid{num}')
        os.makedirs(ogid_outdir, exist_ok=True)
        
        config_paths = {
            0: cells,
            1: [f"{indir}{cell}.table" for cell in cells],
            2: [f"{outdir}ogid{num}/{cell}" for cell in cells],
            3: [f"{indir}{cell}_allregulators.txt" for cell in cells],
            4: [f"{ogid_dir}/{idx}"] * len(cells),
            5: f"{indir}motif_empty.txt" if not args.motifs else [f"{net_dir}{cell}_network.txt" for cell in cells]
        }
        config_df = pd.DataFrame(config_paths)
        config_name = f"testdata_config{'_motifs' if args.motifs else '_noprior'}{num}.txt"
        config_path = os.path.join(indir, config_name)
        config_df.to_csv(config_path, sep='\t', mode='w', index=None, header=False, float_format='%g')
        print(f"Config saved to: {config_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--filelist', help='file list', type=str, default='')
    parser.add_argument('--regfile', help='regulator file', type=str, default='')
    parser.add_argument('--indir', help='directory of input data', type=str, default='')
    parser.add_argument('--outdir', help='output directory', type=str, default='Results/')
    parser.add_argument('--net_dir', help='network directory', type=str, default='')
    parser.add_argument('--splitgene', help='split all genes into X genes per run(0: no split, >0: X)', type=int, default=0)
    parser.add_argument('--addcellsuffix', help='add cell name to the end of gene for each cell', type=int, default=1)
    parser.add_argument('--motifs', help='add motifs as prior or not', type=bool, default=0)
    args = parser.parse_args()
    main(args)
