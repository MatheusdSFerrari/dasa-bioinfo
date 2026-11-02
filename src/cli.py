
import argparse
from pathlib import Path

import pandas as pd
import pysam

from .gene import resolve_genes
from .vcf_processor import process_vcf
from .bam_validator import validate_variant


VERSION = "0.1.0"


def main():
    parser = argparse.ArgumentParser(
        description="Ferramenta de análise de variantes genômicas - Dasa"
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"Dasa Bioinfo {VERSION}"
    )

    parser.add_argument(
        "--vcf",
        required=True,
        help="Caminho do arquivo VCF de entrada"
    )

    parser.add_argument(
        "--gene",
        required=True,
        help="Arquivo TXT contendo os nomes dos genes"
    )

    parser.add_argument(
        "--bam",
        required=False,
        help="Arquivo BAM para validação das variantes (opcional)"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Caminho do arquivo TSV de saída"
    )

    args = parser.parse_args()

    print(f"VCF informado: {args.vcf}", flush=True)
    print(f"Genes informados: {args.gene}", flush=True)
    print(f"Arquivo de saída: {args.output}", flush=True)

    # Obtém as coordenadas dos genes
    regions = resolve_genes(args.gene)

    # Processa o VCF
    print("Processando variantes do VCF...", flush=True)
    df = process_vcf(args.vcf, regions)

    print(f"Variantes identificadas: {len(df)}", flush=True)

    # Validação opcional no BAM
    if args.bam:
        print(f"BAM informado: {args.bam}", flush=True)
        print("Avaliando suporte no BAM...", flush=True)

        results = []

        with pysam.AlignmentFile(args.bam, "rb") as bam:
            if not bam.has_index():
                parser.error(
                    "O arquivo BAM precisa estar indexado."
                )

            for row in df.itertuples(index=False):
                result = validate_variant(
                    bam,
                    row.CHROM,
                    row.POS,
                    row.REF,
                    row.ALT
                )

                results.append(result)

        bam_df = pd.DataFrame(
            results,
            columns=[
                "BAM_RAW_COUNT",
                "BAM_DP",
                "BAM_REF",
                "BAM_ALT",
                "BAM_OTHER",
                "BAM_VAF",
                "BAM_ALT_MAPQ",
                "BAM_ALT_BQ",
                "BAM_STATUS"
            ]
        )

        df = pd.concat(
            [df.reset_index(drop=True), bam_df],
            axis=1
        )

        print("\nResumo da validação:")
        print(df["BAM_STATUS"].value_counts(dropna=False))

    else:
        print(
            "Nenhum BAM informado. "
            "Executando somente a análise do VCF.",
            flush=True
        )

    # Cria o diretório de saída
    output_path = Path(args.output)
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    # Exporta o TSV
    df.to_csv(
        output_path,
        sep="\t",
        index=False
    )

    print(f"Resultado salvo em: {output_path}", flush=True)


if __name__ == "__main__":
    main()
