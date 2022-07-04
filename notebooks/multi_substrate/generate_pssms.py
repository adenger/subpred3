import sys
sys.path.append("../../src")
from main.transporter_dataset import create_dataset
from main.cdhit import cd_hit
from main.eval import (
    optimize_hyperparams,
    preprocess_pandas,
    get_confusion_matrix,
    get_classification_report,
    get_independent_test_set,
    quick_test,
    full_test,
)
from main.compositions import calculate_aac, calculate_paac
from main.pssm import calculate_pssms_notebook
import seaborn as sns

LOG_FILE = "../../logs/athaliana_amino_sugar_potassium.log"


sodium = ["Q2UVJ5", "O04034", "Q9FKP1", "Q9LJI2", "Q9SYG9", "Q84TI7"]
gdu = ["O81775", "Q9SW07", "Q9FHH5", "Q8S8A0", "Q3E965", "Q3EAV6", "Q3E8L0"]
df_all = create_dataset(
    keywords_substrate_filter=[
        "Amino-acid transport",
        "Electron transport",
        "Sugar transport",
        "Potassium transport",
    ],
    keywords_component_filter=["Transmembrane"],
    keywords_transport_filter=["Transport"],
    input_file="../../data/raw/swissprot/uniprot-reviewed_yes.tab.gz",
    multi_substrate="integrate",
    outliers=gdu + sodium,
    verbose=True,
    tax_ids_filter=[3702],
    output_log=LOG_FILE,
    sequence_clustering=70,
)

calculate_pssms_notebook(df_all.sequence, n_threads=78)
