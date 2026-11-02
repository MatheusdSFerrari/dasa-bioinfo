
import gzip
import re
from pathlib import Path


ANNOTATION = (
    Path(__file__).resolve().parent.parent
    / "references"
    / "gencode.v50.annotation.gtf.gz"
)


def resolve_genes(gene_file):
    with open(gene_file, "r", encoding="utf-8") as file:
        genes = [
            line.strip()
            for line in file
            if line.strip() and not line.lstrip().startswith("#")
        ]

    if not genes:
        raise ValueError("A lista de genes está vazia.")

    genes = list(dict.fromkeys(genes))
    matches = {gene: [] for gene in genes}

    with gzip.open(ANNOTATION, "rt") as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9 or fields[2] != "gene":
                continue

            attributes = dict(
                re.findall(r'(\w+)\s+"([^"]+)"', fields[8])
            )

            gene_name = attributes.get("gene_name")

            if gene_name not in matches:
                continue

            matches[gene_name].append({
                "gene": gene_name,
                "chrom": fields[0],
                "start": int(fields[3]) - 1,
                "end": int(fields[4]),
                "ensembl_id": attributes.get("gene_id")
            })

    regions = []

    for gene in genes:
        found = matches[gene]

        if not found:
            raise ValueError(
                f"Gene não encontrado na anotação: {gene}"
            )

        if len(found) != 1:
            raise ValueError(
                f"Símbolo ambíguo na anotação: {gene}"
            )

        region = found[0]
        regions.append(region)

        print(
            f"{gene}: {region['chrom']}:"
            f"{region['start'] + 1}-{region['end']}"
        )

    return regions
