import os
import pandas as pd
import sys
import numpy as np
import argparse
import re

def genOGid(frames,cells,indir,regdf,splitgene):
    ## 1) generate OGid file:
    gene=frames[0]
    for i in range(1,len(cells)):
        gene=gene.merge(frames[i],how='outer',on='Gene')
    gene.fillna('None',inplace=True)
    n=gene.shape[0]
    Gene_OGID='OG'+np.arange(1,n+1).astype(str).astype('object')+'_1'
    NAME=gene.iloc[:,1]
    for i in range(2,len(cells)+1):
        NAME=NAME+','+gene.iloc[:,i]
    ogid={'Gene_OGID':Gene_OGID, 'NAME':NAME}
    ogid = pd.DataFrame(ogid, columns=['Gene_OGID','NAME'])
    outname=indir+'testdata_ogids.txt'
    print('OGID file:',outname)
    ogid.to_csv(outname, index = None, header=True,  sep='\t', mode='w',float_format='%g')
    
    ## 2) TFs_OGs.txt (get OGID for regulators!!)
    regdata=regdf[0]
    for i in range(1,len(cells)):
        regdata=regdata.merge(regdf[i],how='outer',on='Gene')
    regulators=regdata['Gene']
    ## find OG id of regulators (OGID=index+1):
    tfs=gene[gene['Gene'].isin(regulators)==True].index+1
    ## double check if regulator IDs are correct:
    #print("Are regulator IDs correct?",all(gene.loc[tfs-1]['Gene']==regulators))
    tfs=pd.DataFrame({'TF':tfs})
    tffilename=indir+'TFs_OGs.txt'
    print('regulator OG id file:',tffilename)
    tfs.to_csv(tffilename, index = None, header=None,  sep='\t', mode='w',float_format='%g')
    
    ## 3) AllGenes.txt (OGID for all genes)
    filename=indir+'AllGenes.txt'
    print('gene OG id file:',filename)
    geneid=list(range(1,n+1))
    geneid=pd.DataFrame({'gene':geneid})
    geneid.to_csv(filename, index = None, header=None,  sep='\t', mode='w',float_format='%g')

    ## split AllGenes.txt into 50 genes per run
    if splitgene:
        j=0
        outdir=indir+'ogids/'
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        for i in range(1,n+1,splitgene):
            ogs=list(range(i,min(i+splitgene,n+1)))
            #print(ogs)
            ogs=pd.DataFrame({'gene':ogs})
            filename=outdir+'AllGenes'+str(j)+'.txt'
            ogs.to_csv(filename, index = None, header=None,  sep='\t', mode='w',float_format='%g')
            j+=1
        print("Split genes into ogs 0 ->",j)


def main(args):
    # Set up input/output directories
    indir = os.path.join(args.indir, '')
    outdir = os.path.join(args.outdir, '')
    net_dir = os.path.join(args.net_dir, '')
    regfile = args.regfile
    splitgene = args.splitgene
    cellsuffix = args.addcellsuffix
    
    # Read in the files list
    files = pd.read_table(args.filelist, header=None)
    cells = files[0].values  # Cell type names
    filelist = [os.path.join(indir, f'{d}.table') for d in files[0].values]   # Corresponding file paths
    
    
    # 0) Generate celltype_order.txt
    files[0].to_csv(os.path.join(indir, 'celltype_order.txt'), sep='\n', mode='w', index=None, header=False, float_format='%g')
    
    frames = []
    if len(regfile) > 0:
        reg = pd.read_table(regfile, header=None)

    regdf = []
    for i, cell in enumerate(cells):
        infile = filelist[i]
        print(f"Processing {infile} for cell {cell}")

        data = pd.read_table(infile, sep='\t', index_col=0)
        allgenes = data.index.str.split(f'_{cell}').str.join('')

        if cellsuffix == 1:
            d = pd.DataFrame({'Gene': allgenes, cell: allgenes + f'_{cell}'})
        else:
            d = pd.DataFrame({'Gene': allgenes, cell: allgenes})
        
        frames.append(d)
        
        ## 1) Save allGenes.txt
        filename = os.path.join(indir, f'{cell}_allGenes.txt')
        d[cell].to_csv(filename, sep='\n', mode='w', index=None, header=False, float_format='%g')
        print(f"{cell} - Number of genes: {d.shape[0]}") 
        
        ## 2) Save allregulators.txt
        if len(regfile) > 0:
            dreg = reg.loc[reg[0].isin(allgenes)]
        else:
            dreg = pd.DataFrame({0: [gene for gene in allgenes if gene != 'Expression']})
        
        ## Check if all regulators are in the data
        print(f"{cell} - Are filtered regulators all in data? {all(dreg[0].isin(allgenes))}")
        print(f"{cell} - Number of filtered regulators: {dreg.shape[0]}")   
        if cellsuffix == 1:
            df = pd.DataFrame({'Gene': dreg[0].values, cell: dreg[0].values + f'_{cell}'})
        else:
            df = pd.DataFrame({'Gene': dreg[0].values, cell: dreg[0].values})
        regdf.append(df)
        df.to_csv(os.path.join(indir, f"{cell}_allregulators.txt"), sep='\n', mode='w', index=None, header=False, float_format='%g')
        
        ## 3) Update Gene column in data and overwrite file
        if not all(data.index == allgenes + f'_{cell}') and cellsuffix == 1:
            data['Gene'] = allgenes + f'_{cell}'
            data = data.reindex(columns=['Gene'] + list(data.columns[:-1]))
            outname = os.path.join(indir, f"{cell}.table")
            print(f"Overwriting data: {outname}")
            data.to_csv(outname, sep='\t', mode='w', index=None, header=True, float_format='%g')

    # Generate multiple config files based on gene splits
    genOGid(frames, cells, indir, regdf, splitgene)
    
    ogid_dir = os.path.join(indir, 'ogids')
    sub = [f for f in os.listdir(f'{indir}/ogids/') if os.path.isfile(os.path.join(f'{indir}/ogids/', f))]
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
    parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--filelist', help='file list', type=str, default='')
    parser.add_argument('--regfile', help='regulator file', type=str, default='')
    parser.add_argument('--indir', help='directory of input data', type=str, default='')
    parser.add_argument('--outdir', help='output directory', type=str, default='Results/')
    parser.add_argument('--motifs', help='add motifs as prior or not', type=bool, default=0)
    parser.add_argument('--net_dir', help='network directory', type=str, default='')
    parser.add_argument('--splitgene', help='split all genes into X genes per run(0: no split, >0: X)', type=int, default=0)
    parser.add_argument('--addcellsuffix', help='add cell name to the end of gene for each cell: gene_cellname', type=int, default=1)
    args = parser.parse_args()
    main(args)

    parser.add_argument('--motifs', help='add motifs as prior or not', type=bool, default=0)
    args = parser.parse_args()
    main(args)
