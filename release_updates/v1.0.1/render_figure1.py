from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'code'))
import make_expanded_manuscript_figures as ef
plt.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42,'ps.fonttype':42})
maps={'phenotype_manifest.csv':'Figure1B_expanded_phenotypes.csv'}
source=HERE/'figure1_source_data'
def read(path):
    import pandas as pd
    name=Path(path).name
    if name=='phenotype_manifest.csv':
        # Original isolate-level input is included separately, not reconstructed from totals.
        return pd.read_csv(source/name,keep_default_na=False)
    mapping={'assembly_qc_with_checkm2.csv':'Figure1G_expanded_checkm2.csv','panel_callability.csv':'Figure1D_expanded_callability.csv','qc_callability_samples.csv':'Figure1EF_expanded_qc_callability_samples.csv'}
    return pd.read_csv(source/mapping[name],keep_default_na=False)
ef.read_csv=read
ef.write_source=lambda *args: None
def save(fig,folder,stem):
    import numpy as np
    assert len(fig.axes[3].lines)==1
    assert np.all(np.asarray(fig.axes[3].lines[0].get_ydata())==90)
    fig.axes[3].lines[0].remove()
    assert len(fig.axes[6].lines)==2
    fig.axes[3].set_ylabel('Target-callable isolates (%)')
    fig.axes[2].set_ylabel('Materialized genomes')
    for t in fig.axes[0].texts:t.set_text(t.get_text().replace('linked n=','materialized n='))
    fig.axes[6].legend(frameon=False,fontsize=7,loc='lower left')
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/'Figure_1.pdf',bbox_inches='tight')
    fig.savefig(folder/'Figure_1.png',dpi=300,bbox_inches='tight')
    plt.close(fig)
ef.save_figure=save
ef.figure1(HERE,HERE/'rendered_figure1',source)
