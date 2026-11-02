# Desafio Dasa Bioinfo

Ferramenta CLI desenvolvida em Python para identificar variantes de genes de interesse em um arquivo VCF e, opcionalmente, avaliar o suporte aos alelos de referência e alternativo em um arquivo BAM indexado.

A análise é direcionada aos genes FANCA, HEXA e TGFBR2, mas aceita outros símbolos gênicos presentes no arquivo de anotação incluído no projeto.

## Funcionalidades

- Leitura de um VCF e identificação das variantes que se sobrepõem às coordenadas dos genes informados em um arquivo TXT.
- Consulta local à anotação GENCODE v50, baseada no genoma de referência GRCh38, sem necessidade de acesso à API do Ensembl durante a execução.
- Extração de genótipo, zigosidade, profundidade e demais informações disponíveis no VCF.
- Avaliação opcional de evidências no BAM, incluindo contagem de observações REF/ALT, qualidade de mapeamento e qualidade de base.
- Distinção entre variantes sem alinhamentos na região e variantes com alinhamentos, mas sem observações válidas segundo os critérios aplicados.
- Exportação dos resultados para arquivo TSV.
- Execução reprodutível em ambiente Docker.

## Requisitos

- Docker instalado e em funcionamento.
- Arquivo VCF compatível com a referência GRCh38 e com a nomenclatura cromossômica da anotação (`chr1`, `chr2`, etc.).
- Arquivo TXT com os nomes dos genes de interesse separados por linha.
- Opcionalmente, BAM ordenado por coordenadas e respectivo índice `.bai`, com a mesma referência e nomenclatura cromossômica do VCF.

As dependências Python são instaladas na imagem a partir do `requirements.txt`. Entre as principais bibliotecas utilizadas estão pysam e pandas.

## Estrutura do projeto

```text
desafio_bioprod/
├── Dockerfile                  # documento para a reprodução em ambiente Docker
├── requirements.txt            # dependências Python para serem instaladas
├── README.md                   # informações sobre a ferramenta
├── AI_USAGE.md                 # documentação do uso de IA
├── references/                 # Arquivo contendo as informações dos genes de interesse para a extração de coordenadas
│   └── gencode.v50.annotation.gtf.gz
├── src/
│   ├── __init__.py
│   ├── cli.py                 # interface de linha de comando
│   ├── gene.py                # leitura da lista de genes e consulta ao GTF
│   ├── vcf_processor.py       # extração das variantes do VCF
│   └── bam_validator.py       # avaliação opcional das evidências no BAM
└── data/                      # dados locais para a utilização da ferramenta
    ├── input/
    │   ├── genes.txt
    │   ├── desafio.vcf
    │   ├── desafio.bam
    │   └── desafio.bam.bai
    └── output/
```

## Preparação das entradas

Crie as pastas locais de entrada e saída, caso não existam:

```bash
mkdir -p data/input data/output
```

Exemplo de `data/input/genes.txt`:

```text
FANCA
HEXA
TGFBR2
```

O programa também pode receber outros genes presentes no GTF. O arquivo de anotação é lido localmente, e os intervalos são obtidos a partir das entradas do tipo `gene`.

## Construção da imagem Docker

Na raiz do projeto, execute:

```bash
docker build -t dasa-bioinfo:0.1.0 .
```

Para visualizar a versão:

```bash
docker run dasa-bioinfo:0.1.0 --version
```

Para consultar os argumentos disponíveis:

```bash
docker run dasa-bioinfo:0.1.0 --help
```

### Execução 1 — Apenas VCF

```bash
docker run \
  -v "$(pwd)/data:/data" \
  dasa-bioinfo:0.1.0 \
  --vcf /data/input/desafio.vcf \
  --gene /data/input/genes.txt \
  --output /data/output/variantes_vcf.tsv
```

Essa execução identifica variantes nos intervalos gênicos e produz um TSV.

### Execução 2 — VCF com avaliação opcional do BAM

Coloque o BAM e o BAI em `data/input/`, com nomes compatíveis, por exemplo `desafio.bam` e `desafio.bam.bai`.

```bash
docker run \
  -v "$(pwd)/data:/data" \
  dasa-bioinfo:0.1.0 \
  --vcf /data/input/desafio.vcf \
  --gene /data/input/genes.txt \
  --bam /data/input/desafio.bam \
  --output /data/output/variantes_bam.tsv
```

A opção `--bam`; quando informada, acrescenta colunas de avaliação dos alinhamentos ao TSV. O índice do BAM é necessário para acessar regiões genômicas.

### Argumentos da CLI

| Argumento | Obrigatório | Descrição |
|---|---|---|
| `--vcf` | Sim | Caminho do arquivo VCF de entrada. |
| `--gene` | Sim | Caminho do TXT com os nomes dos genes. |
| `--output` | Sim | Caminho do arquivo TSV de saída. |
| `--bam` | Não | Caminho do BAM para avaliação de suporte das variantes. |
| `--version` | Não | Exibe a versão da ferramenta. |
| `--help` | Não | Exibe o menu de ajuda da CLI. |

## Interpretação dos resultados

### Colunas provenientes do VCF

| Coluna | Interpretação |
|---|---|
| `GENE` | Gene cujo intervalo contém ou se sobrepõe à variante. |
| `CHROM` | Cromossomo. |
| `POS` | Posição genômica do registro VCF (1-based). |
| `REF` / `ALT` | Alelos de referência e alternativo. |
| `GT` | Genótipo representado no VCF, incluindo separador de fase quando disponível. |
| `ZYGOSITY` | Classificação de zigosidade derivada do genótipo. |
| `DP` | Profundidade informada no campo da amostra do VCF. |
| `REF_DEPTH` / `ALT_DEPTH` | Contagens de profundidade alélica extraídas de `AD`, quando disponíveis. |
| `VAF` | Fração alélica alternativa calculada a partir das profundidades extraídas de `AD`. |
| `QUAL` | Qualidade da chamada no VCF. |
| `FILTER` | Filtros indicados no registro VCF. |

Atenção: na implementação atual, `ALT_DEPTH` e `VAF` são calculados a partir do primeiro alelo alternativo. Para registros multialélicos, a interpretação desses campos requer cautela.

### Colunas adicionais quando `--bam` é informado

| Coluna | Interpretação |
|---|---|
| `BAM_RAW_COUNT` | Número de registros de alinhamento que se sobrepõem ao intervalo consultado, antes dos filtros de validação. |
| `BAM_DP` | Total de observações classificadas como REF, ALT ou OTHER após os filtros aplicados. |
| `BAM_REF` | Observações que sustentam o alelo de referência. |
| `BAM_ALT` | Observações que sustentam o alelo alternativo. |
| `BAM_OTHER` | Observações válidas que não foram classificadas como REF nem ALT. |
| `BAM_VAF` | Fração `BAM_ALT / BAM_DP`. |
| `BAM_ALT_MAPQ` | Média da qualidade de mapeamento das observações ALT aceitas. |
| `BAM_ALT_BQ` | Média da qualidade de base das observações ALT aceitas; para indels, sua interpretação depende da base âncora e das bases verificadas pelo algoritmo. |
| `BAM_STATUS` | Resultado da avaliação de evidências na região. |

A versão padrão utiliza os limiares MAPQ ≥ 20 e BQ ≥ 20 e exclui, entre outros, alinhamentos secundários, suplementares, duplicados e marcados como falha de controle de qualidade. O processamento de pares sobrepostos pelo `pysam.pileup()` pode reduzir a contagem em comparação com os alinhamentos brutos.

#### Valores possíveis de `BAM_STATUS`

| Status | Significado |
|---|---|
| `SUPORTE_ALT_OBSERVADO` | Há pelo menos uma observação do alelo ALT aceita pelos filtros. |
| `ALT_NAO_OBSERVADO` | Há observações válidas, mas nenhuma sustenta ALT. |
| `SEM_ALINHAMENTOS` | Não há alinhamentos sobrepostos ao intervalo consultado no BAM fornecido. |
| `SEM_COBERTURA_VALIDA` | Há alinhamentos, mas nenhuma observação foi aceita pelos filtros. |
| `NAO_AVALIAVEL` | O tipo ou a representação da variante não é contemplado pelos filtros de avaliação. |

Campos vazios nas colunas BAM indicam que a métrica não foi calculada.

### Limitações importantes

- Interpretação clínica: a presença de suporte ALT não constitui, isoladamente, confirmação clínica, classificação de patogenicidade ou diagnóstico.
- Compatibilidade das referências: diferenças de assembly e de nomenclatura cromossômica podem levar à ausência de resultados ou comparações incorretas.

## Reprodutibilidade e integridade dos dados

- Utilize a mesma versão da imagem Docker, do `requirements.txt` e da anotação `gencode.v50.annotation.gtf.gz` para reproduzir uma execução.
- Preserve os arquivos originais de entrada e seus checksums.
- O processamento de genes utiliza uma anotação local, sem consultas externas durante a execução.

## Referências técnicas

- [pysam — documentação](https://pysam.readthedocs.io/en/latest/api.html)
- [GENCODE — anotações do genoma humano](https://www.gencodegenes.org/human/)
- [SAM/BAM/VCF — especificações de formatos (hts-specs)](https://github.com/samtools/hts-specs)
