O ChatGPT (OpenAI) foi utilizado como ferramenta de apoio à implementação do software desenvolvido para o desafio técnico. As principais consultas envolveram a anotação GENCODE (GRCh38), além da implementação de funções para avaliação automatizada do suporte aos alelos em arquivos BAM, considerando proporções de leituras e critérios de qualidade de mapeamento (MAPQ) e de base (BQ).

A principal estratégia de validação consistiu na execução individual dos comandos sugeridos pela IA diretamente no terminal Linux, antes de sua incorporação ao código definitivo. Os resultados foram comparados com as informações dos arquivos de entrada e com o comportamento esperado de cada operação.

Quando identificados erros ou inconsistências, foram realizadas consultas e comparações a repositórios e documentações técnicas, seguidas de ajustes e novas execuções. Após a verificação individual, as funcionalidades foram integradas aos módulos Python e à interface de linha de comando (CLI), sendo posteriormente avaliadas em conjunto.

Durante o desenvolvimento, a instabilidade da API REST do Ensembl demonstrou a necessidade de adotar uma alternativa, que foi sugerida pela IA, baseada na consulta local à anotação GENCODE (GRCh38), posteriormente testada e incorporada à ferramenta.

A validação inicial dos alelos foi realizada por meio da inspeção visual dos alinhamentos no Integrative Genomics Viewer (IGV), considerando as proporções de leituras que sustentavam os alelos de referência e alternativo, além das métricas de qualidade de mapeamento (MAPQ) e de base (BQ). Para automatizar essa avaliação e integrá-la à ferramenta, foi utilizado o auxílio da IA na implementação de funções em Python capazes de aplicar os critérios de qualidade e quantificar o suporte aos alelos diretamente no arquivo BAM. Os resultados obtidos foram posteriormente comparados com as evidências observadas no IGV.

As sugestões da IA foram submetidas à avaliação crítica e a testes manuais antes de sua utilização. A interpretação dos resultados e as decisões finais de implementação permaneceram sob responsabilidade do desenvolvedor.
