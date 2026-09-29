import pandas as pd 
from pathlib import Path
import matplotlib.pyplot as plt
BASE_DIR=Path(__file__).resolve().parent.parent
DATA_DIR=BASE_DIR/"data"
FIGURES_DIR = BASE_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True) 
df = pd.read_csv(
    DATA_DIR / "brca1_clinvar.tsv",
    sep="\t",
    index_col=False,
    usecols=range(24)
)
print(df.head())
print(df.columns)
print(df.shape)
with open(DATA_DIR / "brca1_clinvar.tsv", "r", encoding="utf-8") as file:
    first_line = file.readline()
    second_line = file.readline()

print(df.shape)
print(df.index[:5])
print(df["Germline classification"].head())
print(df["Germline review status"].head())
print(df["Germline classification"].value_counts(dropna=False,))
classification_counts=(df["Germline classification"].value_counts(dropna=False)) 
classification_counts.plot(kind="bar")
plt.title("BRCA1 Germline classififcation Distribution")
plt.xlabel("Germline classififcation")
plt.ylabel("Number of Variants")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(
    FIGURES_DIR / "germline_classification_distribution.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()
variant_type_counts=df["Variant type"].value_counts(dropna=False)
print(variant_type_counts)
variant_type_counts.plot(kind="bar")
plt.title("BRCA1 Variant type Distribution")
plt.xlabel("Variant type")
plt.ylabel("Number of variant")
plt.xticks(rotation=45,ha="right")
plt.tight_layout()
plt.savefig(
    FIGURES_DIR / "variant_type_distribution.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()
molecular_counts=df["Molecular consequence"].value_counts(dropna=False)
print(molecular_counts)
molecular_split=(
    df["Molecular consequence"]
    .dropna()
    .str.split("|")
    .explode()
)
print(molecular_split.value_counts())
top_molecular=molecular_split.value_counts().head(10)
top_molecular.plot(kind="bar")
plt.title("Top molecular consequences in BRCA1")
plt.xlabel("Molecular consequence")
plt.ylabel("Count")
plt.xticks(rotation=45,ha="right")
plt.tight_layout()
plt.savefig(
    FIGURES_DIR / "molecular_consequence_distribution.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()
pathogenic_variants=df[
    df["Germline classification"]=="Pathogenic"
]
pathogenic_by_type=pathogenic_variants["Variant type"].value_counts()
total_by_type=df["Variant type"].value_counts()
pathogenic_pourcentage_by_type=pathogenic_by_type/ total_by_type*100
print(pathogenic_pourcentage_by_type.sort_values(ascending=False))
summary_by_type=pd.DataFrame({
    "Total":total_by_type,
    "Pathogenic":pathogenic_by_type
})
summary_by_type["Pathogenic"]=summary_by_type["Pathogenic"].fillna(0)
summary_by_type["Pathogenic%"]=(summary_by_type["Pathogenic"]/summary_by_type["Total"]*100)
print(summary_by_type.sort_values("Pathogenic%",ascending=False))
reliable_type=summary_by_type[summary_by_type["Total"]>=100]
print(
    reliable_type.sort_values(
        "Pathogenic%",
        ascending=False
     
    )
)
reliable_type["Pathogenic%"].plot(kind="bar")
plt.title("Pathogenic Percentage by BRCA1 Variant Type")
plt.xlabel("Variant type")
plt.ylabel("percentage")
plt.xticks(rotation=45,ha="right")
plt.tight_layout()
plt.savefig(
    FIGURES_DIR / "pathogenic_percentage_by_variant_type.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()
pathogenic_example=df[
    (df["Germline classification"]=="Pathogenic")
    & (df["Molecular consequence"].notna())
    & (df["Protein change"].notna())
].iloc[0]
print(pathogenic_example[
    [
        "Name",
        "Protein change",
        "Variant type",
        "Molecular consequence",
        "Germline classification",
        "Germline review status"
    ]
    
    
])
pathogenic_example=df[
    (df["Germline classification"]=="Pathogenic")
    & (df["Molecular consequence"].notna())
    & (df["Protein change"].notna())
].iloc[0]
print(pathogenic_example[
    [
        "Name",
        "Protein change",
        "Variant type",
        "Molecular consequence",
        "Germline classification",
        "Germline review status"
    ]
    
    
])
vus_examples = df[
    (df["Germline classification"] == "Uncertain significance")
    & (df["Molecular consequence"].notna())
    & (df["Protein change"].notna())
]

vus_example = vus_examples.iloc[0]

print(vus_example[
    [
        "Name",
        "Protein change",
        "Variant type",
        "Molecular consequence",
        "Germline classification",
        "Germline review status"
    ]
])
benign_examples = df[
    (df["Germline classification"] == "Benign")
    & (df["Molecular consequence"].notna())
    & (df["Protein change"].notna())
]

benign_example = benign_examples.iloc[0]

print(benign_example[
    [
        "Name",
        "Protein change",
        "Variant type",
        "Molecular consequence",
        "Germline classification",
        "Germline review status"
    ]
])

print("\n--- Pathogenic Example ---")
print(pathogenic_example[
    ["Name", "Variant type", "Molecular consequence", "Germline review status"]
])

print("\n--- VUS Example ---")
print(vus_example[
    ["Name", "Variant type", "Molecular consequence", "Germline review status"]
])
print("\n--- Dataset Summary ---")
print("Total variants:", len(df))
print("Pathogenic variants:", (df["Germline classification"] == "Pathogenic").sum())
print("VUS variants:", (df["Germline classification"] == "Uncertain significance").sum())
print("Benign variants:", (df["Germline classification"] == "Benign").sum())
print("\n--- Benign Example ---")
print(benign_example[
    ["Name", "Variant type", "Molecular consequence", "Germline review status"]
])

reliable_type.to_csv(
    BASE_DIR / "reliable_variant_type_summary.csv",
    index=True
)

print("\nAnalysis completed successfully.")
print("Results saved in:", BASE_DIR)

