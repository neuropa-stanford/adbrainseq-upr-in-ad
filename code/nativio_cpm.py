# CPM-normalized Nativio log2FC (fixes library-size global shift)
import math, os, openpyxl
BASE=("/data/adbrainseq/Publication/2024_scRNA seq/"
      "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication")
def load_nativio_cpm():
    wb=openpyxl.load_workbook(os.path.join(BASE,"SuppD3_NativioSeq.xlsx"),read_only=True); ws=wb.active
    rows=[]; genes=[]
    for r in ws.iter_rows(min_row=3,values_only=True):
        g=r[1]
        if not g or str(g)=='None': continue
        try: ad=[float(x) for x in r[6:18]]; old=[float(x) for x in r[18:28]]
        except: continue
        genes.append(str(g).strip()); rows.append(ad+old)
    wb.close()
    import numpy as np
    M=np.array(rows); tot=M.sum(axis=0)             # per-sample library size
    CPM=M/tot*1e6
    ad=CPM[:,:12].mean(axis=1); old=CPM[:,12:].mean(axis=1)
    out={}
    for i,g in enumerate(genes):
        if ad[i]>0 and old[i]>0:
            lf=math.log2(ad[i]/old[i])
            if g not in out or abs(lf)>0: out[g]=lf
    return out
if __name__=='__main__':
    import statistics
    d=load_nativio_cpm()
    print('Nativio CPM: %d genes, global median log2FC = %+.4f'%(len(d),statistics.median(d.values())))
