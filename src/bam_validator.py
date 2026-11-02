
import pysam


MIN_MAPQ = 20
MIN_BQ = 20


def validate_variant(bam, chrom, pos, ref, alt):
    ref_count = 0
    alt_count = 0
    other_count = 0
    alt_mapq = []
    alt_bq = []

    # Nesta versão, não avaliamos variantes complexas.
    if not ref or not alt or "," in alt:
        return empty_result("NAO_AVALIAVEL")

    # Identifica o tipo de variante
    if len(ref) == 1 and len(alt) == 1:
        variant_type = "SNV"
    elif len(ref) > 1 and len(alt) == 1 and ref[0] == alt:
        variant_type = "DEL"
    elif len(ref) == 1 and len(alt) > 1 and alt[0] == ref:
        variant_type = "INS"
    else:
        return empty_result("NAO_AVALIAVEL")

    pos0 = int(pos) - 1
    ref = ref.upper()
    alt = alt.upper()

    # Verifica se existem alinhamentos brutos na posição
    raw_count = bam.count(
        chrom,
        pos0,
        pos0 + len(ref),
        read_callback="nofilter"
    )

    if raw_count == 0:
        result = empty_result("SEM_ALINHAMENTOS")
        result["BAM_RAW_COUNT"] = 0
        return result

    for column in bam.pileup(
        chrom,
        pos0,
        pos0 + 1,
        truncate=True,
        stepper="all",
        min_base_quality=0,
        min_mapping_quality=0,
        ignore_overlaps=True,
        max_depth=100000
    ):
        if column.reference_pos != pos0:
            continue

        for pileup_read in column.pileups:
            read = pileup_read.alignment

            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.is_duplicate
                or read.is_qcfail
                or read.mapping_quality < MIN_MAPQ
                or pileup_read.is_del
                or pileup_read.is_refskip
            ):
                continue

            qpos = pileup_read.query_position
            sequence = read.query_sequence
            qualities = read.query_qualities

            if (
                qpos is None
                or sequence is None
                or qualities is None
                or qualities[qpos] < MIN_BQ
            ):
                continue

            allele = "OTHER"
            bq = qualities[qpos]

            if variant_type == "SNV":
                base = sequence[qpos].upper()

                if base == ref:
                    allele = "REF"
                elif base == alt:
                    allele = "ALT"

            elif variant_type == "DEL":
                deletion_size = len(ref) - len(alt)

                if (
                    sequence[qpos].upper() == ref[0]
                    and pileup_read.indel == -deletion_size
                ):
                    allele = "ALT"

                elif pileup_read.indel == 0:
                    # Verifica se a leitura contém
                    # toda a sequência de referência.
                    pairs = read.get_aligned_pairs(
                        matches_only=True
                    )

                    positions = {
                        rpos: q
                        for q, rpos in pairs
                    }

                    bases = []
                    base_qualities = []

                    for rpos in range(pos0, pos0 + len(ref)):
                        q = positions.get(rpos)

                        if q is None:
                            break

                        bases.append(sequence[q].upper())
                        base_qualities.append(qualities[q])

                    if (
                        len(bases) == len(ref)
                        and "".join(bases) == ref
                        and min(base_qualities) >= MIN_BQ
                    ):
                        allele = "REF"

            elif variant_type == "INS":
                insertion_size = len(alt) - len(ref)

                if (
                    sequence[qpos].upper() == ref
                    and pileup_read.indel == insertion_size
                ):
                    inserted = sequence[
                        qpos + 1:qpos + 1 + insertion_size
                    ].upper()

                    inserted_qualities = qualities[
                        qpos:qpos + 1 + insertion_size
                    ]

                    if (
                        inserted == alt[1:]
                        and len(inserted_qualities) == insertion_size + 1
                        and min(inserted_qualities) >= MIN_BQ
                    ):
                        allele = "ALT"
                        bq = min(inserted_qualities)

                elif pileup_read.indel == 0:
                    # Exige que a leitura continue após
                    # a base âncora da inserção.
                    positions = read.get_reference_positions(
                        full_length=False
                    )

                    if (
                        sequence[qpos].upper() == ref
                        and pos0 + 1 in positions
                    ):
                        allele = "REF"

            if allele == "REF":
                ref_count += 1

            elif allele == "ALT":
                alt_count += 1
                alt_mapq.append(read.mapping_quality)
                alt_bq.append(bq)

            else:
                other_count += 1

    depth = ref_count + alt_count + other_count

    if depth == 0:
        status = "SEM_COBERTURA_VALIDA"
    elif alt_count > 0:
        status = "SUPORTE_ALT_OBSERVADO"
    else:
        status = "ALT_NAO_OBSERVADO"

    return {
	"BAM_RAW_COUNT": raw_count,
        "BAM_DP": depth,
        "BAM_REF": ref_count,
        "BAM_ALT": alt_count,
        "BAM_OTHER": other_count,
        "BAM_VAF": round(alt_count / depth, 4) if depth else None,
        "BAM_ALT_MAPQ": (
            round(sum(alt_mapq) / len(alt_mapq), 2)
            if alt_mapq else None
        ),
        "BAM_ALT_BQ": (
            round(sum(alt_bq) / len(alt_bq), 2)
            if alt_bq else None
        ),
        "BAM_STATUS": status
    }


def empty_result(status):
    return {
        "BAM_RAW_COUNT": None,
	"BAM_DP": None,
        "BAM_REF": None,
        "BAM_ALT": None,
        "BAM_OTHER": None,
        "BAM_VAF": None,
        "BAM_ALT_MAPQ": None,
        "BAM_ALT_BQ": None,
        "BAM_STATUS": status
    }
