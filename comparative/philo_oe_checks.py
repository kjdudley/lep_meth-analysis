#!/usr/bin/env python3
"""The three checks on Philopotamus's anomalous CpG o/e (DNMT_MOTIF_SCALE.md).

Check 1 (per-gene distribution) was done locally and came back UNIMODAL, just
shifted, and matched-GC strata did not explain it either. So the two live
hypotheses are:

  (a) genuinely less CpG depletion, i.e. recent methylation in Annulipalpia;
  (b) gene SPANS inflated by CpG-undepleted repeat in introns, which would
      shift every gene uniformly and therefore also look unimodal.

EXONS-ONLY is what separates them: exons carry little repeat. If exon-only o/e
converges on the other species, it was (b). If it stays high, (a) survives.

Soft-mask state is reported because if the assembly carries lowercase repeat
annotation the masked recompute is free; if not, that check needs a real
repeat library and is deferred rather than faked.

o/e is (nCpG / (nC * nG)) * n, identical to cpg_oe.py:67 so the numbers are
comparable to the gene-span run.
"""
import sys, gzip, collections

def opener(p): return gzip.open(p,"rt") if p.endswith(".gz") else open(p)

def load_fa(p):
    seqs, name, buf = {}, None, []
    with opener(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name: seqs[name]="".join(buf)
                name=line[1:].split()[0]; buf=[]
            else: buf.append(line.strip())
    if name: seqs[name]="".join(buf)
    return seqs

def counts(s, unmasked_only=False):
    # CHECK 2. Soft-masked (lowercase) bases are replaced by N rather than
    # deleted: deleting them would splice distant uppercase bases together and
    # manufacture CpGs that do not exist. N counts as neither C nor G and
    # cannot form a CpG, which is what we want.
    u = "".join(c if c.isupper() else "N" for c in s) if unmasked_only else s.upper()
    return u.count("C"), u.count("G"), u.count("CG"), len(u)-u.count("N")

def oe(nC,nG,nCpG,n):
    return (nCpG/(nC*nG))*n if nC and nG and n else None

def main(label, fa_path, gff_path):
    seqs=load_fa(fa_path)
    soft=sum(1 for s in seqs.values() for c in s[:2000000] if c.islower())
    sampled=sum(min(len(s),2000000) for s in seqs.values())
    print(f"  soft-masked fraction (first 2 Mb/seq): {soft/sampled:.4f}"
          f"  {'-> repeat-masked recompute possible' if soft/sampled>0.02 else '-> NOT soft-masked, masked check deferred'}")

    genes=collections.defaultdict(list); exons=collections.defaultdict(list)
    with opener(gff_path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            f=line.rstrip("\n").split("\t")
            if len(f)<9: continue
            if f[2]=="gene":  genes[f[0]].append((int(f[3])-1,int(f[4])))
            elif f[2]=="exon": exons[f[0]].append((int(f[3])-1,int(f[4])))

    for what, d, unmasked in (("gene spans",genes,False),("gene spans UNMASKED",genes,True),
                              ("EXONS only",exons,False),("EXONS UNMASKED",exons,True)):
        C=G=CG=N=0; nfeat=0
        # merge overlaps so shared exons are not double counted
        for sq,iv in d.items():
            s=seqs.get(sq)
            if not s: continue
            iv.sort(); merged=[]
            for a,b in iv:
                if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
                else: merged.append((a,b))
            for a,b in merged:
                c,g,cg,n=counts(s[a:b], unmasked); C+=c; G+=g; CG+=cg; N+=n; nfeat+=1
        v=oe(C,G,CG,N)
        print(f"  {what:<12} merged features={nfeat:<7} bp={N/1e6:8.1f} Mb  "
              f"GC={(C+G)/N:.4f}  CpG o/e = {v:.5f}" if v else f"  {what}: undefined")

if __name__=="__main__":
    main(*sys.argv[1:4])
