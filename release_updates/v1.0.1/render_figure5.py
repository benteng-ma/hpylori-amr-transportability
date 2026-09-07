"""Re-render Figure 5 from unchanged repository result tables. No model fitting."""
from pathlib import Path
import sys, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'code'))
import make_expanded_manuscript_figures as ef
parser=argparse.ArgumentParser()
parser.add_argument('--repository-root',type=Path,default=HERE.parents[1])
args=parser.parse_args()
ROOT=args.repository_root.resolve()
UP=HERE/"rendered_figure5"
UP.mkdir(exist_ok=True)
OUT=HERE
WORK=UP

def fig5():
    plt.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42,'ps.fonttype':42})
    def save(fig, output_dir, stem):
        folds=ef.numeric(ef.read_csv(ROOT/'results/lineage_validation/leakage_benchmark_folds.csv'), ['balanced_accuracy'])
        models=['MUTATION_ONLY_LOGISTIC','LINEAGE_ONLY_LOGISTIC','MUTATION_PLUS_LINEAGE_LOGISTIC']
        labels=['Mutation','Lineage','Mutation + lineage']
        colors=['#2878B5','#E39C37','#7A5195']
        groups=[('RANDOM_ISOLATE_SPLIT',None,'Random'),('CLONE_GROUPED_SPLIT',None,'Clone'),('LEAVE_COUNTRY_OUT',None,'Country'),('LEAVE_LINEAGE_OUT',None,'Lineage'),('HPGP_TO_EXTERNAL','CHINA_NINGXIA_2022','Ningxia'),('HPGP_TO_EXTERNAL','ZENODO_10369064','Reads')]
        rng=np.random.default_rng(20260831)
        for ax,drug,letter in zip(fig.axes[:2],['clarithromycin','levofloxacin'],['A','B']):
            ax.clear()
            for mi,model in enumerate(models):
                for si,(split,cohort,_) in enumerate(groups):
                    sub=folds[(folds.antibiotic==drug)&(folds.model==model)&(folds.split_type==split)]
                    if cohort:
                        sub=sub[sub.fold==cohort]
                    vals=sub.balanced_accuracy.dropna().to_numpy()
                    if not len(vals): continue
                    x=si+(mi-1)*.22
                    if len(vals)>1:
                        ax.scatter(x+rng.normal(0,.022,len(vals)),vals,s=11,alpha=.28,color=colors[mi])
                    ax.errorbar(x,vals.mean(),yerr=vals.std(ddof=1) if len(vals)>1 else 0,fmt='D',color=colors[mi],capsize=3,markersize=4)
            ax.axvline(3.5,color='#C6C6C6',lw=.7)
            ax.axhline(.5,color='#AAAAAA',ls=':')
            ax.set_xticks(np.arange(6),[g[2] for g in groups],rotation=30,ha='right')
            ax.set_ylim(0,1.04); ax.set_xlim(-.5,5.5); ax.set_ylabel('Balanced accuracy')
            if drug=='clarithromycin':
                ax.text(4,.08,'NE',ha='center',fontsize=8,color='#666666')
            ef.style_axis(ax); ef.panel_label(ax,letter)
        fig.axes[0].legend(handles=[Line2D([0],[0],marker='D',ls='',color=c,label=l) for c,l in zip(colors,labels)],frameon=False,fontsize=7,ncol=3,loc='lower left',bbox_to_anchor=(0,1.015),borderaxespad=0,columnspacing=.9,handletextpad=.35)
        fig.savefig(UP/'Figure_5.pdf',bbox_inches='tight')
        fig.savefig(UP/'Figure_5.png',dpi=300,bbox_inches='tight')
        plt.close(fig)
        ext=folds[folds.split_type=='HPGP_TO_EXTERNAL']
        ext.to_csv(WORK/'Figure5AB_external_cohort_values.csv',index=False)
    ef.save_figure=save
    old_write=ef.write_source
    ef.write_source=lambda frame, dest, name: old_write(frame,dest,name.replace('Figure4','Figure5'))
    ef.figure4(ROOT,UP,WORK/'figure5_source_data')

fig5()
