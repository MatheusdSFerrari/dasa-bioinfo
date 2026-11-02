
import pysam
import pandas as pd


def process_vcf(vcf_path, regions):
    results = []

    with pysam.VariantFile(vcf_path) as vcf:
        for variant in vcf:
            for region in regions:

                if not (
                    variant.chrom == region["chrom"]
                    and variant.start < region["end"]
                    and variant.stop > region["start"]
                ):
                    continue

                sample = next(iter(variant.samples.values()))

                gt = sample.get("GT")
                dp = sample.get("DP")
                ad = sample.get("AD")

                # Formatação do genótipo
                if gt is None:
                    genotype = "./."
                    zygosity = "Sem chamada"
                else:
                    separator = "|" if sample.phased else "/"
                    genotype = separator.join(
                        "." if allele is None else str(allele)
                        for allele in gt
                    )

                    called = [a for a in gt if a is not None]

                    if len(called) != len(gt):
                        zygosity = "Chamada parcial"
                    elif len(set(called)) > 1:
                        zygosity = "Heterozigoto"
                    elif called[0] == 0:
                        zygosity = "Homozigoto referência"
                    else:
                        zygosity = "Homozigoto alternativo"

                # Profundidade dos alelos
                ref_depth = ad[0] if ad and len(ad) > 0 else None
                alt_depth = ad[1] if ad and len(ad) > 1 else None

                # Frequência do alelo alternativo
                if (
                    ref_depth is not None
                    and alt_depth is not None
                    and ref_depth + alt_depth > 0
                ):
                    vaf = alt_depth / (ref_depth + alt_depth)
                else:
                    vaf = None

                results.append({
                    "GENE": region["gene"],
                    "CHROM": variant.chrom,
                    "POS": variant.pos,
                    "REF": variant.ref,
                    "ALT": ",".join(variant.alts or []),
                    "GT": genotype,
                    "ZYGOSITY": zygosity,
                    "DP": dp,
                    "REF_DEPTH": ref_depth,
                    "ALT_DEPTH": alt_depth,
                    "VAF": round(vaf, 4) if vaf is not None else None,
                    "QUAL": round(variant.qual, 2)
                            if variant.qual is not None else None,
                    "FILTER": ";".join(variant.filter.keys())
                })

    return pd.DataFrame(results)
