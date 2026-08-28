import glob, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

TRUTH_DIR=os.path.expanduser("~/foldback/sim_data")
RESULTS_DIR=os.path.expanduser("~/foldback/results")
PROBE_DIR=os.path.expanduser("~/foldback/results/read_level_detector/probe")
METHODS=["yacrd","breakinator","read_level_detector_probe"]
FOLDS=["quarter","middle","near_end"]
FRACTIONS=[1,5,10]

def load_calls(method,condition):
    if method=="read_level_detector_probe":
        df=pd.read_csv(f"{PROBE_DIR}/calls_read_level_detector_{condition}.tsv",sep="\t")
        df["flagged"]=pd.to_numeric(df["flagged"])>0.8
    else:
        path=glob.glob(f"{RESULTS_DIR}/**/calls_{method}_{condition}.tsv",recursive=True)[0]
        df=pd.read_csv(path,sep="\t")
        df["flagged"]=df["flagged"].astype(str).str.strip().eq("True")
    return df.groupby("read_id",as_index=False)["flagged"].max()

def load_truth(condition):
    t=pd.read_csv(f"{TRUTH_DIR}/truth_{condition}.tsv",sep="\t")
    t["is_foldback"]=t["is_foldback"].astype(str).str.strip().eq("True")
    if "adapter_present" in t:
        t["adapter_present"]=t["adapter_present"].astype(str).str.strip().eq("True")
    return t

def merged(method,condition):
    d=load_truth(condition).merge(load_calls(method,condition),on="read_id",how="left")
    d["flagged"]=d["flagged"].fillna(False).astype(bool)
    return d

# TP/TN/FP/FN for every simulated condition
rows=[]
for method in METHODS:
    for fold in FOLDS:
        for frac in FRACTIONS:
            condition=f"{frac}pct_{fold}"
            d=merged(method,condition)
            tp=(d.is_foldback&d.flagged).sum()
            tn=(~d.is_foldback&~d.flagged).sum()
            fp=(~d.is_foldback&d.flagged).sum()
            fn=(d.is_foldback&~d.flagged).sum()
            rows.append([method,fold,frac,tp,tn,fp,fn,tp/(tp+fn) if tp+fn else np.nan])

cm=pd.DataFrame(rows,columns=["method","fold","fraction","TP","TN","FP","FN","recall"])
print("\n=== SIMULATED CONDITIONS ===")
print(cm.to_string(index=False))

# Clean-control confusion matrix + FPR
rows=[]
for method in METHODS:
    d=merged(method,"clean_control")
    tp=(d.is_foldback&d.flagged).sum()
    tn=(~d.is_foldback&~d.flagged).sum()
    fp=(~d.is_foldback&d.flagged).sum()
    fn=(d.is_foldback&~d.flagged).sum()
    fpr=fp/(fp+tn) if fp+tn else np.nan
    rows.append([method,tp,tn,fp,fn,fpr])

clean=pd.DataFrame(rows,columns=["method","TP","TN","FP","FN","FPR"])
print("\n=== CLEAN CONTROL ===")
print(clean.to_string(index=False))

# Adapter-specific recall
rows=[]
for method in METHODS:
    for fold in FOLDS:
        for frac in FRACTIONS:
            d=merged(method,f"{frac}pct_{fold}")
            for adapter in [True,False]:
                pos=d.is_foldback&(d.adapter_present==adapter)
                tp=(pos&d.flagged).sum()
                fn=(pos&~d.flagged).sum()
                rows.append([method,fold,frac,adapter,tp,fn,tp/(tp+fn) if tp+fn else np.nan])

df=pd.DataFrame(rows,columns=["method","fold","fraction","adapter_present","TP","FN","recall"])

FOLDS=["middle","quarter","near_end"]
FOLD_LABELS={"middle":"Middle","quarter":"Off-center","near_end":"Near-end"}
fpr=dict(zip(clean["method"],clean["FPR"]))

fig,axes=plt.subplots(3,3,figsize=(18,11),sharey=True)
x=np.arange(len(METHODS)); w=.42

for r,fold in enumerate(FOLDS):
    for c,frac in enumerate(FRACTIONS):
        ax=axes[r,c]
        s=df[(df.fold==fold)&(df.fraction==frac)]
        present=[s[(s.method==m)&s.adapter_present].recall.iloc[0]*100 for m in METHODS]
        absent=[s[(s.method==m)&~s.adapter_present].recall.iloc[0]*100 for m in METHODS]

        ax.bar(x-w/2,present,w,label="Adapter present")
        ax.bar(x+w/2,absent,w,label="Adapter absent")

        for i,m in enumerate(METHODS):
            y=max(present[i],absent[i])+3
            ax.text(i,y,f"FPR {fpr[m]*100:.1f}%",ha="center",fontsize=12)

        ax.set_title(f"{frac}% - {FOLD_LABELS[fold]}",loc="left",fontweight="bold",fontsize=16)
        ax.set_ylim(0,115)
        ax.set_xticks(x)
        ax.set_xticklabels(METHODS,fontsize=12)
        ax.tick_params(axis="y",labelsize=12)
        ax.spines[["top","right"]].set_visible(False)

        if c==0:
            ax.set_ylabel("Recall (%)",fontsize=14)

h,l=axes[0,0].get_legend_handles_labels()
fig.legend(h,l,loc="upper left",bbox_to_anchor=(.04,.96),ncol=2,frameon=False,fontsize=14)
fig.suptitle("Framing A",fontsize=26,fontweight="bold",y=1.01)
fig.tight_layout(rect=[0,0,1,.95])
plt.show()
