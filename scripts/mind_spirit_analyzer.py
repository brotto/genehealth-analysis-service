"""
Analisador de Mente & Espírito
Analisa variantes genéticas associadas a traços de personalidade, predisposições
de saúde mental e sensibilidade espiritual (incluindo genética da mediunidade).

TRÊS CATEGORIAS:
  1. Traços de Personalidade - COMT, CADM2, MSRA, OPRM1, OXTR
  2. Saúde Mental & Resposta ao Estresse - FKBP5, BDNF, MTHFR, CACNA1C, RGS2, SNAP25
  3. Sensibilidade Espiritual - HTR2A, VMAT2/SLC18A2, OXTR (contexto de absorção)

BASE CIENTÍFICA:
  - Personalidade: SNPs validados por GWAS (Nature Human Behaviour 2024, Translational Psychiatry 2023)
  - Saúde Mental: Estudos do consórcio PGC (Cell 2024, Nature Genetics 2025)
  - Espiritual: Estudo USP/BJPsych 2025 sobre mediunidade (PMID 39874024), pesquisa HTR2A sobre experiências místicas

RESSALVAS IMPORTANTES:
  - Todas as associações são probabilísticas, não determinísticas
  - Os tamanhos de efeito são geralmente pequenos (cada SNP explica <1% da variância do traço)
  - Ambiente, cultura e escolha pessoal modificam fortemente qualquer tendência genética
  - Relatórios de saúde mental são apenas educativos, não diagnósticos
  - Achados sobre sensibilidade espiritual são preliminares e exploratórios

References:
  - PGC MDD Working Group 2024, Cell (PMC11092713) - 697 depression loci
  - Anxiety GWAS 2025, Nature Genetics - 58 anxiety loci
  - Wagner et al. 2025, BJPsych (PMID 39874024) - mediumship candidate genes
  - Pasternak & Zangari 2025, Questao de Ciencia - critical review
  - Singh et al. 2023, Translational Psychiatry - CADM2 impulsivity
"""

from typing import Dict, List, Tuple, Any, Optional


# ---------------------------------------------------------------------------
# BANCO DE SNPs: Traços de Personalidade
# ---------------------------------------------------------------------------

PERSONALITY_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs4680",
        "gene": "COMT",
        "name": "Regulação de Dopamina e Estilo Cognitivo",
        "chromosome": "22",
        "position": 19963748,
        "effectSize": "moderate",
        "genotypes": {
            "GG": {
                "label": "Val/Val",
                "nickname": "Guerreiro",
                "interpretation": (
                    "Maior atividade da enzima COMT leva a uma eliminação mais rápida de dopamina no "
                    "córtex pré-frontal. Associado a melhor tolerância ao estresse, maior limiar de dor "
                    "e vantagem sob pressão. Pode apresentar desempenho cognitivo basal ligeiramente menor, "
                    "mas melhor desempenho sob estresse."
                ),
                "score": 0.7,
            },
            "AG": {
                "label": "Val/Met",
                "nickname": "Equilibrado",
                "interpretation": (
                    "Atividade intermediária da COMT proporciona um perfil de dopamina equilibrado. "
                    "Você provavelmente tem um estilo cognitivo versátil — tolerância adequada ao estresse "
                    "com boa flexibilidade cognitiva basal. Genótipo mais comum."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "Met/Met",
                "nickname": "Preocupado",
                "interpretation": (
                    "Menor atividade da enzima COMT significa eliminação mais lenta de dopamina, resultando "
                    "em níveis mais altos de dopamina pré-frontal. Associado a desempenho cognitivo superior "
                    "no estado basal (melhor memória de trabalho, atenção), mas potencialmente maior "
                    "vulnerabilidade ao estresse e ansiedade sob pressão."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "COMT Val158Met (rs4680) é um dos SNPs funcionais mais estudados na genética comportamental. "
            "O alelo Val (G) produz uma enzima de alta atividade que elimina dopamina ~4x mais rápido "
            "que o alelo Met (A). Meta-análises confirmam efeitos pequenos mas reais na resposta ao "
            "estresse (d ≈ 0,15-0,25) e na flexibilidade cognitiva. O modelo 'Guerreiro vs Preocupado' "
            "(Goldman et al. 2005) captura o trade-off: Val/Val se destaca sob estresse mas tem menor "
            "dopamina pré-frontal basal, enquanto Met/Met se destaca na cognição basal mas é mais "
            "sensível ao estresse."
        ),
        "references": [
            {"pmid": "16151010", "title": "The 'warrior' and 'worrier' model for COMT", "year": 2005},
            {"pmid": "29520078", "title": "COMT Val158Met and cognition meta-analysis", "year": 2018},
        ],
        "actionableInsights": [
            "Val/Val: Você pode prosperar em ambientes de alta pressão — considere aproveitar isso nas escolhas de carreira",
            "Met/Met: Práticas de mindfulness e redução de estresse podem ser especialmente benéficas para você",
            "Exercício físico regular ajuda a otimizar os níveis de dopamina independentemente do genótipo",
        ],
    },
    {
        "rsid": "rs17518584",
        "gene": "CADM2",
        "name": "Tolerância ao Risco e Impulsividade",
        "chromosome": "3",
        "position": 85890189,
        "effectSize": "moderate",
        "genotypes": {
            "CC": {
                "label": "Típico",
                "nickname": "Cauteloso",
                "interpretation": (
                    "Genótipo de referência associado a níveis típicos de tolerância ao risco. "
                    "Você provavelmente tem uma abordagem equilibrada para decisões de risco."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "Intermediário",
                "nickname": "Moderado",
                "interpretation": (
                    "Uma cópia do alelo associado ao risco. Tendência ligeiramente elevada "
                    "para comportamento de risco e busca por novidades em comparação com portadores CC."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "Elevado",
                "nickname": "Tomador de risco",
                "interpretation": (
                    "Duas cópias do alelo associado ao risco. Associado a maior tolerância "
                    "ao risco, busca por novidades e impulsividade em estudos GWAS. "
                    "Também relacionado a tendências de uso de substâncias em estudos populacionais."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "CADM2 (Cell Adhesion Molecule 2) rs17518584 é um dos achados GWAS mais bem validados "
            "para comportamento de risco. Replicado no UK Biobank (n>400.000) e nas coortes de pesquisa "
            "da 23andMe. CADM2 é expresso no cérebro e desempenha um papel no desenvolvimento sináptico. "
            "Os tamanhos de efeito são pequenos mas consistentes (Singh et al. 2023, Translational "
            "Psychiatry). Um estudo de 2023 da Nature Human Behaviour identificou 27+ loci independentes "
            "para tolerância ao risco, com CADM2 entre os mais fortes."
        ),
        "references": [
            {"pmid": "36658148", "title": "CADM2 implicated in impulsive personality by GWAS", "year": 2023},
            {"pmid": "30643258", "title": "Risk tolerance GWAS in UK Biobank", "year": 2019},
        ],
        "actionableInsights": [
            "Portadores TT: Estar ciente das tendências elevadas de impulsividade pode ajudar em decisões financeiras e de saúde",
            "Considere usar estruturas de tomada de decisão para escolhas importantes da vida",
            "Canalize tendências de risco de forma construtiva através de empreendedorismo ou esportes",
        ],
    },
    {
        "rsid": "rs4925638",
        "gene": "MSRA",
        "name": "Irritabilidade e Inibição Comportamental",
        "chromosome": "8",
        "position": 10040702,
        "effectSize": "moderate",
        "genotypes": {
            "AA": {
                "label": "Típico",
                "nickname": "Temperamento estável",
                "interpretation": (
                    "Genótipo de referência associado a níveis típicos de reatividade emocional "
                    "e inibição comportamental."
                ),
                "score": 0.3,
            },
            "AG": {
                "label": "Intermediário",
                "nickname": "Reatividade moderada",
                "interpretation": (
                    "Uma cópia do alelo variante. Reatividade emocional ligeiramente elevada "
                    "em comparação com portadores AA."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "Elevado",
                "nickname": "Reatividade elevada",
                "interpretation": (
                    "Duas cópias do alelo variante. Associado a maior irritabilidade "
                    "e redução da inibição comportamental em estudos GWAS. MSRA desempenha um papel "
                    "na proteção contra estresse oxidativo no cérebro."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "MSRA (Methionine Sulfoxide Reductase A) rs4925638 foi associado a irritabilidade, "
            "tolerância ao risco e impulsividade em análises GWAS do UK Biobank. MSRA protege proteínas "
            "contra danos oxidativos — variantes podem afetar a resiliência neuronal ao estresse e vias "
            "de regulação emocional. Replicado em Davies et al. 2017 (PMC5537199)."
        ),
        "references": [
            {"pmid": "28696412", "title": "Replication of CADM2 and MSRA on human behavior", "year": 2017},
        ],
        "actionableInsights": [
            "Portadores GG: Dieta rica em antioxidantes (frutas vermelhas, folhas verdes) pode apoiar a função da MSRA",
            "Exercício físico regular é um dos mais fortes moderadores da irritabilidade",
            "Técnicas cognitivo-comportamentais podem ajudar a gerenciar a reatividade emocional elevada",
        ],
    },
    {
        "rsid": "rs1799971",
        "gene": "OPRM1",
        "name": "Sensibilidade Emocional e Processamento de Recompensa",
        "chromosome": "6",
        "position": 154360797,
        "effectSize": "moderate",
        "genotypes": {
            "AA": {
                "label": "Asp/Asp",
                "nickname": "Sensibilidade típica",
                "interpretation": (
                    "Genótipo de referência para o receptor mu-opioide. Sensibilidade à dor padrão "
                    "e resposta placebo típica. Genótipo mais comum mundialmente."
                ),
                "score": 0.4,
            },
            "AG": {
                "label": "Asp/Asn",
                "nickname": "Sensibilidade aumentada",
                "interpretation": (
                    "Uma cópia do alelo G (A118G). Associado a sensibilidade alterada à dor, "
                    "resposta placebo reduzida e potencialmente maior sensibilidade à frustração. "
                    "Pode vivenciar rejeição social de forma mais intensa."
                ),
                "score": 0.6,
            },
            "GG": {
                "label": "Asn/Asn",
                "nickname": "Alta sensibilidade",
                "interpretation": (
                    "Duas cópias do alelo variante. Maior sensibilidade emocional a interações "
                    "sociais e processamento de recompensa potencialmente alterado. Menos responsivo "
                    "a efeitos placebo e pode necessitar de abordagens diferentes para manejo da dor."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "OPRM1 A118G (rs1799971) afeta o receptor mu-opioide, que medeia dor, recompensa e "
            "vínculo social. O alelo G reduz a expressão do receptor e foi associado a maior "
            "sensibilidade à frustração, resposta placebo reduzida e processamento alterado de "
            "dor social. Múltiplos estudos de neuroimagem mostram que o alelo G está associado a "
            "maiores respostas neurais à rejeição social (Way et al. 2009, PNAS)."
        ),
        "references": [
            {"pmid": "19622738", "title": "OPRM1 A118G and social rejection sensitivity", "year": 2009},
            {"pmid": "15583718", "title": "Functional effects of the OPRM1 A118G polymorphism", "year": 2005},
        ],
        "actionableInsights": [
            "Portadores G: Sua sensibilidade emocional elevada pode ser uma força em relacionamentos empáticos",
            "Esteja ciente de que a percepção de dor e a resposta a medicamentos podem diferir das médias populacionais",
            "Conexão social e redes de apoio podem ser especialmente importantes para portadores G",
        ],
    },
    {
        "rsid": "rs53576",
        "gene": "OXTR",
        "name": "Sensibilidade Social e Empatia",
        "chromosome": "3",
        "position": 8804371,
        "effectSize": "weak",
        "genotypes": {
            "GG": {
                "label": "GG",
                "nickname": "Maior tendência empática",
                "interpretation": (
                    "Associado a maiores escores de empatia, maior sensibilidade social "
                    "e reatividade emocional mais forte a sinais sociais em estudos iniciais. "
                    "Nota: replicações em larga escala mostraram resultados mistos — os tamanhos "
                    "de efeito são muito pequenos."
                ),
                "score": 0.7,
            },
            "AG": {
                "label": "AG",
                "nickname": "Intermediário",
                "interpretation": (
                    "Sensibilidade intermediária do receptor de ocitocina. Tendências típicas "
                    "de vínculo social. Genótipo mais comum."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "AA",
                "nickname": "Menor sensibilidade social",
                "interpretation": (
                    "Associado a sensibilidade social e escores de empatia ligeiramente menores "
                    "em alguns estudos. No entanto, replicações em larga escala no UK Biobank "
                    "(n>100.000) encontraram efeitos fracos e inconsistentes."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "OXTR rs53576 é a variante do receptor de ocitocina mais estudada para comportamento social. "
            "Estudos iniciais associaram o genótipo GG a maior empatia e comportamento pró-social. "
            "No entanto, uma meta-análise de 2017 (n>10.000) encontrou a associação fraca e "
            "inconsistente (PMID 28343138). Análises do UK Biobank mostram efeito mínimo. Este SNP "
            "é incluído por completude, mas os resultados devem ser interpretados com cautela."
        ),
        "references": [
            {"pmid": "28343138", "title": "Meta-analysis of OXTR rs53576 and empathy", "year": 2017},
            {"pmid": "21775986", "title": "OXTR and social behavior review", "year": 2011},
        ],
        "actionableInsights": [
            "O comportamento social é predominantemente moldado pelo ambiente e experiências pessoais",
            "Empatia é uma habilidade que pode ser desenvolvida independentemente do genótipo",
            "Conexões sociais fortes beneficiam todos — invista em relacionamentos",
        ],
    },
    {
        "rsid": "rs2254298",
        "gene": "OXTR",
        "name": "Vínculo Social e Apego",
        "chromosome": "3",
        "position": 8798609,
        "effectSize": "weak",
        "genotypes": {
            "GG": {
                "label": "GG",
                "nickname": "Referência",
                "interpretation": (
                    "Genótipo de referência para esta variante do OXTR. Padrões típicos de "
                    "apego e vínculo social."
                ),
                "score": 0.4,
            },
            "GA": {
                "label": "GA",
                "nickname": "Intermediário",
                "interpretation": (
                    "Uma cópia do alelo A. Alguns estudos associam a padrões alterados de "
                    "vínculo social, mas os achados são inconsistentes entre populações."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "AA",
                "nickname": "Variante",
                "interpretation": (
                    "Duas cópias do alelo A. Associado a estilos de apego alterados "
                    "em alguns estudos, particularmente em populações do leste asiático. Replicação "
                    "limitada em coortes de ascendência europeia."
                ),
                "score": 0.6,
            },
        },
        "scientificBasis": (
            "OXTR rs2254298 foi estudado para associações com vínculo social, segurança de apego "
            "e comportamento pró-social. Os resultados variam significativamente por população — o "
            "alelo A mostra efeitos mais fortes em estudos com populações do leste asiático do que em "
            "amostras de ascendência europeia. Valor preditivo limitado para comportamento social individual."
        ),
        "references": [
            {"pmid": "21775986", "title": "OXTR polymorphisms and social behavior", "year": 2011},
        ],
        "actionableInsights": [
            "Os estilos de apego são formados principalmente pelas experiências de cuidado na infância",
            "O apego seguro pode ser desenvolvido em qualquer idade através de terapia e relacionamentos",
        ],
    },
]


# ---------------------------------------------------------------------------
# BANCO DE SNPs: Saúde Mental e Resposta ao Estresse
# ---------------------------------------------------------------------------

MENTAL_HEALTH_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs1360780",
        "gene": "FKBP5",
        "name": "Resposta ao Estresse e Regulação do Eixo HPA",
        "chromosome": "6",
        "position": 35639794,
        "effectSize": "strong",
        "genotypes": {
            "CC": {
                "label": "CC",
                "nickname": "Resposta ao estresse padrão",
                "interpretation": (
                    "Genótipo de referência. Expressão normal de FKBP5 e feedback do receptor "
                    "de glicocorticoides. Resposta e recuperação padrão do cortisol ao estresse."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "CT",
                "nickname": "Sensibilidade ao estresse moderadamente aumentada",
                "interpretation": (
                    "Uma cópia do alelo T de risco. Expressão moderadamente aumentada de FKBP5 "
                    "após estresse, o que pode desacelerar a recuperação do feedback de cortisol. "
                    "Sensibilidade levemente elevada ao estresse, particularmente no contexto de "
                    "experiências adversas."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "TT",
                "nickname": "Sensibilidade ao estresse aumentada",
                "interpretation": (
                    "Duas cópias do alelo T (genótipo de alta indução). Aumento significativo "
                    "na expressão do mRNA de FKBP5 após exposição ao estresse, prejudicando "
                    "o feedback do receptor de glicocorticoides e prolongando a resposta do cortisol. "
                    "Associação mais rigorosa com TEPT e depressão relacionada ao estresse "
                    "no contexto de adversidade na infância."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "FKBP5 rs1360780 é um dos achados mais robustos de SNP único na psiquiatria relacionada "
            "ao estresse. O alelo T cria uma resposta de metilação de DNA alelo-específica ao estresse "
            "que aumenta a expressão da proteína FKBP5, que se liga e inibe o receptor de glicocorticoides, "
            "perturbando o feedback do cortisol. Uma revisão de 2025 da MDPI Genes (Zannas et al.) "
            "caracterizou fenótipos comportamentais. A interação gene×ambiente com maus-tratos na "
            "infância é bem replicada em múltiplas coortes independentes. Os tamanhos de efeito são "
            "modestos mas mecanisticamente convincentes."
        ),
        "references": [
            {"pmid": "39856776", "title": "FKBP5 rs1360780: genetic variation and behavioral phenotypes", "year": 2025},
            {"pmid": "24029109", "title": "FKBP5 epigenetics and stress vulnerability", "year": 2013},
        ],
        "actionableInsights": [
            "Portadores TT: Práticas de manejo do estresse (meditação, exercício) podem ser especialmente importantes",
            "Intervenção precoce para sintomas relacionados ao estresse é valiosa para portadores do alelo T",
            "Abordagens terapêuticas que visam a reatividade ao estresse (TCC, EMDR) são opções baseadas em evidências",
            "Esta variante NÃO determina sua saúde mental — ambiente e estratégias de enfrentamento importam enormemente",
        ],
    },
    {
        "rsid": "rs6265",
        "gene": "BDNF",
        "name": "Neuroplasticidade e Regulação do Humor",
        "chromosome": "11",
        "position": 27658369,
        "effectSize": "moderate",
        "genotypes": {
            "CC": {
                "label": "Val/Val",
                "nickname": "Neuroplasticidade padrão",
                "interpretation": (
                    "Genótipo de referência (Val66). Secreção normal de BDNF e liberação "
                    "dependente de atividade. Neuroplasticidade e função hipocampal padrão."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "Val/Met",
                "nickname": "Moderadamente alterado",
                "interpretation": (
                    "Uma cópia do alelo Met. Secreção reduzida de BDNF dependente de atividade. "
                    "Efeitos sutis na memória e volume hipocampal em alguns estudos. A associação "
                    "com depressão per se é fraca no nível populacional."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "Met/Met",
                "nickname": "Secreção reduzida de BDNF",
                "interpretation": (
                    "Duas cópias do alelo Met. Redução significativa na liberação de BDNF "
                    "dependente de atividade. Associado a menor volume hipocampal, consolidação "
                    "alterada de memória e, em alguns estudos, risco modestamente elevado de depressão — "
                    "particularmente em homens e populações idosas."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "BDNF Val66Met (rs6265) afeta a secreção dependente de atividade do fator neurotrófico "
            "derivado do cérebro, crucial para a plasticidade sináptica. Uma meta-análise de 2023 da "
            "Frontiers in Psychiatry encontrou efeito quase nulo para TDM (OR ~0,96). No entanto, o "
            "alelo Met mostra associações mais fortes em homens (OR 1,27) e depressão geriátrica "
            "(OR ~1,48). O maior valor está na neuroplasticidade e biologia do exercício — o BDNF "
            "medeia muitos benefícios cognitivos relacionados ao exercício."
        ),
        "references": [
            {"pmid": "37398592", "title": "BDNF Val66Met and depression meta-analysis", "year": 2023},
            {"pmid": "15457404", "title": "BDNF Val66Met and hippocampal function", "year": 2004},
        ],
        "actionableInsights": [
            "Exercício é a forma mais conhecida de aumentar os níveis de BDNF — especialmente importante para portadores Met",
            "Exercício aeróbico (corrida, natação, ciclismo) tem mais evidências para elevação de BDNF",
            "Portadores Met podem se beneficiar particularmente de treinamento cognitivo e experiências de aprendizado novas",
        ],
    },
    {
        "rsid": "rs1801133",
        "gene": "MTHFR",
        "name": "Metabolismo do Folato e Humor",
        "chromosome": "1",
        "position": 11856378,
        "effectSize": "moderate",
        "genotypes": {
            "GG": {
                "label": "CC (normal)",
                "nickname": "Atividade enzimática plena",
                "interpretation": (
                    "Atividade normal da enzima MTHFR (~100%). Metabolismo eficiente do folato "
                    "e ciclo de metilação. Sem preocupações de humor relacionadas ao folato."
                ),
                "score": 0.2,
            },
            "AG": {
                "label": "CT (heterozigoto)",
                "nickname": "Atividade reduzida (~65%)",
                "interpretation": (
                    "Uma cópia do alelo T (C677T). Atividade da enzima MTHFR reduzida para "
                    "aproximadamente 65% do normal. Geralmente clinicamente insignificante, mas "
                    "a ingestão adequada de folato é aconselhável."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "TT (homozigoto)",
                "nickname": "Significativamente reduzida (~30%)",
                "interpretation": (
                    "Duas cópias do alelo T. Atividade da enzima MTHFR reduzida para aproximadamente "
                    "30% do normal. Meta-análise de 30 estudos encontra OR = 1,20 para depressão "
                    "(TT vs CC). Efeito concentrado em populações asiáticas. A prevalência de TT no "
                    "Brasil é de aproximadamente 10%. O valor acionável é alto: suplementação de "
                    "metilfolato para portadores TT é amplamente aceita clinicamente."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "MTHFR C677T (rs1801133) reduz a atividade da metilenotetrahidrofolato redutase, "
            "afetando o ciclo de metabolismo de um carbono que produz SAMe (o principal doador de "
            "metil para síntese de neurotransmissores). Homozigotos TT têm ~30% de atividade "
            "enzimática e podem desenvolver hiperhomocisteinemia leve. Meta-análise mostra associação "
            "modesta com depressão (OR = 1,20, p=0,0004) concentrada em populações asiáticas. O maior "
            "valor clínico é a acionabilidade: suplementação de folato/metilfolato."
        ),
        "references": [
            {"pmid": "23212058", "title": "MTHFR C677T and depression meta-analysis", "year": 2013},
            {"pmid": "15166018", "title": "MTHFR and psychiatric disorders review", "year": 2004},
        ],
        "actionableInsights": [
            "Portadores TT: Considere suplementação de L-metilfolato (consulte um profissional de saúde)",
            "Garanta ingestão adequada de folato na dieta: folhas verdes, leguminosas, alimentos fortificados",
            "Portadores TT devem monitorar os níveis de homocisteína com seu médico",
            "B12 e B6 são cofatores importantes — garanta ingestão adequada",
        ],
    },
    {
        "rsid": "rs1006737",
        "gene": "CACNA1C",
        "name": "Canal de Cálcio e Estabilidade do Humor",
        "chromosome": "12",
        "position": 2345295,
        "effectSize": "strong",
        "genotypes": {
            "GG": {
                "label": "GG",
                "nickname": "Referência",
                "interpretation": (
                    "Genótipo de referência. Função padrão do canal de cálcio tipo L e "
                    "excitabilidade neuronal."
                ),
                "score": 0.2,
            },
            "AG": {
                "label": "AG",
                "nickname": "Um alelo de risco",
                "interpretation": (
                    "Uma cópia do alelo A de risco. Expressão ligeiramente alterada do canal de "
                    "cálcio. Risco transdiagnóstico modestamente elevado para condições de humor "
                    "(bipolar, depressão, esquizofrenia) — OR aproximadamente 1,07."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "AA",
                "nickname": "Dois alelos de risco",
                "interpretation": (
                    "Duas cópias do alelo A. A variante de risco genético mais consistentemente "
                    "replicada para transtorno bipolar (OR ~1,14). Esta é uma variante transdiagnóstica "
                    "— confere pequenos aumentos de risco em múltiplas condições de humor e "
                    "psicóticas. A maioria dos portadores nunca desenvolve qualquer condição."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "CACNA1C rs1006737 é o SNP de risco para transtorno bipolar mais consistentemente replicado "
            "em todas as meta-análises do PGC (OR ~1,14, p = 9,78×10⁻¹⁰). Codifica a subunidade alfa-1C "
            "do canal de cálcio tipo L dependente de voltagem, crítico para excitabilidade neuronal "
            "e plasticidade sináptica. Importante: esta é uma variante transdiagnóstica — também aumenta "
            "o risco para esquizofrenia e depressão recorrente. O estudo PGC4 trans-ancestral de 2025 "
            "encontrou 93 loci bipolares (23 novos), com CACNA1C permanecendo entre os mais fortes."
        ),
        "references": [
            {"pmid": "34002096", "title": "PGC Bipolar GWAS identifies 64 loci", "year": 2021},
            {"pmid": "21926972", "title": "CACNA1C and psychiatric phenotypes", "year": 2011},
        ],
        "actionableInsights": [
            "Portadores AA: Monitoramento do humor pode ajudar a identificar padrões precocemente",
            "Horários regulares de sono são particularmente importantes para a estabilidade do humor",
            "Esta variante aumenta o risco ligeiramente — a maioria dos portadores nunca desenvolve transtornos de humor",
            "Se você tem histórico familiar de transtorno bipolar, converse com seu médico",
        ],
    },
    {
        "rsid": "rs4606",
        "gene": "RGS2",
        "name": "Sensibilidade à Ansiedade e Reatividade ao Estresse",
        "chromosome": "1",
        "position": 192839582,
        "effectSize": "moderate",
        "genotypes": {
            "CC": {
                "label": "CC",
                "nickname": "Referência",
                "interpretation": (
                    "Genótipo de referência. Expressão padrão de RGS2 e regulação da sinalização "
                    "de proteína G. Níveis típicos de resposta à ansiedade."
                ),
                "score": 0.3,
            },
            "CG": {
                "label": "CG",
                "nickname": "Intermediário",
                "interpretation": (
                    "Uma cópia do alelo G. Expressão potencialmente reduzida de RGS2, o que "
                    "pode aumentar sutilmente a sinalização de hormônios do estresse a jusante "
                    "dos receptores de CRH."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "GG",
                "nickname": "Sensibilidade à ansiedade elevada",
                "interpretation": (
                    "Duas cópias do alelo G. Associado a maior sensibilidade à ansiedade "
                    "e reatividade ao estresse em múltiplos estudos. RGS2 amortece a sinalização "
                    "de proteína G — RGS2 reduzido pode levar a resposta aumentada de hormônios "
                    "do estresse."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "RGS2 (Regulator of G-protein Signaling 2) rs4606 modula a sinalização a jusante "
            "dos receptores de hormônio liberador de corticotrofina (CRH). O alelo G foi associado "
            "a ansiedade e reatividade ao estresse em múltiplos estudos independentes (2008-2024). "
            "Os tamanhos de efeito são pequenos mas o achado foi replicado. Um GWAS de ansiedade "
            "marco de 2025 da Nature Genetics (n=122.341) identificou 58 loci, apoiando ainda mais "
            "o envolvimento de vias GABAérgicas e de sinalização de estresse."
        ),
        "references": [
            {"pmid": "18266781", "title": "RGS2 and anxiety disorders", "year": 2008},
        ],
        "actionableInsights": [
            "Portadores GG: Técnicas estruturadas de relaxamento (relaxamento muscular progressivo, respiração profunda) podem ser especialmente úteis",
            "Exercício físico regular é uma das intervenções ansiolíticas mais fortes",
            "Considere se estratégias de manejo da ansiedade beneficiariam sua rotina diária",
        ],
    },
    {
        "rsid": "rs3746544",
        "gene": "SNAP25",
        "name": "Sinalização Sináptica e Atenção",
        "chromosome": "20",
        "position": 10284862,
        "effectSize": "weak",
        "genotypes": {
            "TT": {
                "label": "TT",
                "nickname": "Referência",
                "interpretation": (
                    "Genótipo de referência. Função padrão da proteína SNAP25 e liberação "
                    "de vesículas sinápticas."
                ),
                "score": 0.3,
            },
            "TG": {
                "label": "TG",
                "nickname": "Intermediário",
                "interpretation": (
                    "Uma cópia do alelo G. Efeitos sutis na eficiência da transmissão sináptica. "
                    "Alguns estudos de associação com TDAH mostram elevação modesta de risco."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "GG",
                "nickname": "Variante",
                "interpretation": (
                    "Duas cópias do alelo G. Associado a tendências de TDAH em alguns "
                    "estudos, embora os tamanhos de efeito sejam pequenos. SNAP25 está envolvido na "
                    "exocitose de vesículas sinápticas e liberação de neurotransmissores."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "SNAP25 (Synaptosomal-Associated Protein 25kDa) rs3746544 é um gene candidato "
            "para TDAH. SNAP25 é essencial para a fusão de vesículas sinápticas e liberação de "
            "neurotransmissores. Os tamanhos de efeito são pequenos. Os principais genes candidatos "
            "clássicos para TDAH (DRD4, DAT1 VNTRs) não estão disponíveis em chips de DNA do "
            "consumidor. O GWAS de TDAH do PGC (2023, PMC10914347) identificou 27 loci com "
            "significância genômica ampla através de abordagens GWAS."
        ),
        "references": [
            {"pmid": "15742474", "title": "SNAP25 and ADHD association", "year": 2005},
        ],
        "actionableInsights": [
            "Atenção e foco são habilidades altamente treináveis independentemente do genótipo",
            "Rotinas estruturadas e organização do ambiente podem apoiar a atenção",
            "Se você suspeita de dificuldades de atenção, avaliação profissional é recomendada",
        ],
    },
]


# ---------------------------------------------------------------------------
# BANCO DE SNPs: Sensibilidade Espiritual e Mediunidade
# ---------------------------------------------------------------------------

SPIRITUAL_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs6313",
        "gene": "HTR2A",
        "name": "Receptor de Serotonina 2A e Intensidade de Experiências Místicas",
        "chromosome": "13",
        "position": 47471478,
        "effectSize": "moderate",
        "genotypes": {
            "CC": {
                "label": "C/C (102C)",
                "nickname": "Maior densidade de receptores",
                "interpretation": (
                    "Genótipo de referência associado a maior densidade de receptores 5-HT2A no "
                    "córtex. Intensidade padrão de experiências transcendentes ou místicas. "
                    "Este é o genótipo mais comum."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "C/T (heterozigoto)",
                "nickname": "Densidade intermediária de receptores",
                "interpretation": (
                    "Uma cópia do alelo T (T102C). Densidade intermediária de receptores 5-HT2A. "
                    "Sensibilidade potencialmente moderada a estados transcendentes ou alterados de "
                    "consciência. Pesquisas com psilocibina mostram intensidade intermediária de "
                    "experiências místicas com este genótipo."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "T/T (102T)",
                "nickname": "Menor densidade de receptores",
                "interpretation": (
                    "Duas cópias do alelo T. Associado a menor densidade de receptores 5-HT2A "
                    "no córtex. Paradoxalmente, estudos mostram que indivíduos com menos receptores "
                    "5-HT2A relatam experiências místicas MAIS intensas durante a administração de "
                    "psilocibina. Isso sugere um mecanismo de sensibilidade compensatória — menos "
                    "receptores, mas maior sensibilidade individual de cada receptor."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "HTR2A (receptor de serotonina 5-HT2A) é o alvo principal dos psicodélicos clássicos "
            "e desempenha papel central na consciência, percepção e experiências místicas. "
            "Um estudo de 2021 (PMID 33501857) descobriu que a ligação cerebral do receptor 5-HT2A "
            "prediz a intensidade dos efeitos místicos durante a administração de psilocibina — "
            "pessoas com menos receptores (genótipo TT) paradoxalmente relatam experiências MAIS "
            "intensas. rs6313 (T102C) é uma variante sinônima em forte desequilíbrio de ligação "
            "com a variante promotora rs6311 (A-1438G). O receptor 5-HT2A também está implicado na "
            "absorção — o traço psicológico mais consistentemente associado a experiências espirituais "
            "e mediúnicas relatadas."
        ),
        "references": [
            {"pmid": "33501857", "title": "Brain 5-HT2A receptor binding predicts mystical experiences", "year": 2021},
            {"pmid": "17606772", "title": "HTR2A polymorphisms and hallucinogen response", "year": 2007},
        ],
        "actionableInsights": [
            "Portadores TT: Você pode ter sensibilidade aumentada a práticas contemplativas (meditação, oração)",
            "Práticas de mindfulness e meditação podem produzir experiências particularmente vívidas para portadores TT",
            "Esta variante se relaciona com consciência e percepção — não com qualquer habilidade espiritual específica",
        ],
    },
    {
        "rsid": "rs6311",
        "gene": "HTR2A",
        "name": "Promotor de Serotonina 2A e Experiências Transcendentes",
        "chromosome": "13",
        "position": 47471867,
        "effectSize": "moderate",
        "genotypes": {
            "GG": {
                "label": "G/G (A-1438)",
                "nickname": "Promotor de referência",
                "interpretation": (
                    "Genótipo de referência para o promotor do HTR2A. Níveis padrão de "
                    "transcrição e expressão do receptor 5-HT2A."
                ),
                "score": 0.3,
            },
            "GA": {
                "label": "G/A (heterozigoto)",
                "nickname": "Expressão intermediária",
                "interpretation": (
                    "Uma cópia do alelo A. Atividade alterada do promotor afetando a expressão "
                    "de 5-HT2A. Fenótipo intermediário para sensibilidade serotoninérgica."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "A/A (-1438G>A)",
                "nickname": "Expressão alterada",
                "interpretation": (
                    "Duas cópias do alelo A. Esta variante do promotor altera a transcrição "
                    "de HTR2A. Em desequilíbrio de ligação com rs6313 — juntos, modulam o papel "
                    "do sistema serotoninérgico na consciência e percepção. Associado a resposta "
                    "alterada a estímulos serotoninérgicos."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "HTR2A rs6311 (A-1438G) é um polimorfismo do promotor em forte desequilíbrio de ligação "
            "com rs6313 (T102C). Juntos, modulam a densidade e função do receptor 5-HT2A. O receptor "
            "de serotonina 2A é o principal mediador de experiências psicodélicas e místicas, e "
            "variações neste sistema receptor foram associadas a diferenças individuais em absorção, "
            "abertura a experiências e autotranscendência. Estudos com gêmeos mostram 40-50% de "
            "herdabilidade para autotranscendência (Bouchard et al.)."
        ),
        "references": [
            {"pmid": "11779785", "title": "HTR2A promoter polymorphism and receptor density", "year": 2002},
        ],
        "actionableInsights": [
            "Combinadas com rs6313, essas variantes moldam o perfil de sensibilidade do seu sistema serotoninérgico",
            "Tradições contemplativas de diversas culturas descrevem fenômenos consistentes com variação serotoninérgica",
        ],
    },
    {
        "rsid": "rs4570625",
        "gene": "SLC18A2",
        "name": "Transporte Vesicular de Monoaminas e Autotranscendência",
        "chromosome": "10",
        "position": 119003566,
        "effectSize": "preliminary",
        "genotypes": {
            "GG": {
                "label": "GG",
                "nickname": "Referência",
                "interpretation": (
                    "Genótipo de referência para a região VMAT2. Função padrão do transportador "
                    "vesicular de monoaminas e empacotamento de neurotransmissores."
                ),
                "score": 0.3,
            },
            "GT": {
                "label": "GT",
                "nickname": "Intermediário",
                "interpretation": (
                    "Uma cópia do alelo T. Na hipótese do 'Gene de Deus' de Dean Hamer (2004), "
                    "o alelo T foi associado a maiores escores de autotranscendência no inventário "
                    "de personalidade TCI de Cloninger. Este achado NÃO foi replicado em GWAS "
                    "de larga escala. Interprete com extrema cautela."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "TT",
                "nickname": "Variante",
                "interpretation": (
                    "Duas cópias do alelo T. O estudo original de Hamer (n~1.000) sugeriu "
                    "maior autotranscendência. No entanto, a hipótese do 'Gene de Deus' NÃO "
                    "foi replicada em estudos de larga escala. VMAT2 é um transportador real de "
                    "monoaminas envolvido no empacotamento de dopamina, serotonina e norepinefrina "
                    "em vesículas — mas vinculá-lo especificamente à espiritualidade permanece "
                    "cientificamente não validado."
                ),
                "score": 0.6,
            },
        },
        "scientificBasis": (
            "VMAT2 (SLC18A2) foi proposto como o 'Gene de Deus' por Dean Hamer em 2004. A variante "
            "rs4570625 foi associada a escores de autotranscendência em uma amostra pequena "
            "(n~1.000). VMAT2 empacota monoaminas (dopamina, serotonina, norepinefrina) em vesículas "
            "sinápticas — tem um papel genuíno na neurotransmissão. No entanto, a hipótese do 'Gene "
            "de Deus' NÃO foi replicada em GWAS de larga escala. O consenso científico considera isso "
            "especulativo. Incluímos por completude e interesse cultural, claramente marcado como "
            "evidência preliminar."
        ),
        "references": [
            {"pmid": "15520385", "title": "Hamer: The God Gene (book review)", "year": 2004},
        ],
        "actionableInsights": [
            "Autotranscendência é uma dimensão de personalidade validada independentemente de sua base genética",
            "Práticas contemplativas (meditação, oração, imersão na natureza) cultivam a transcendência em todas as pessoas",
            "O 'Gene de Deus' é um conceito cultural, não um achado científico estabelecido",
        ],
    },
]


# ---------------------------------------------------------------------------
# CONTEXTO DO ESTUDO DA USP SOBRE MEDIUNIDADE
# ---------------------------------------------------------------------------

USP_STUDY_CONTEXT = (
    "Em janeiro de 2025, pesquisadores afiliados à USP (Universidade de São Paulo) e ao "
    "IPq (Instituto de Psiquiatria) publicaram um estudo marco no Brazilian Journal "
    "of Psychiatry (PMID 39874024): 'Candidate genes related to spiritual mediumship: "
    "a whole-exome sequencing analysis of highly gifted mediums.' O estudo comparou o "
    "sequenciamento completo do exoma de 54 médiuns altamente experientes (10+ anos de prática) "
    "contra 53 parentes de primeiro grau não-médiuns. Principais achados: 15.669 variantes "
    "genéticas foram encontradas exclusivamente em médiuns, com 33 genes alterados em ≥1/3 "
    "dos médiuns mas em nenhum de seus parentes. A via biológica mais afetada foi a do sistema "
    "inflamatório e imune (43,9%), com a translocação de ZAP-70 para a sinapse imunológica "
    "como o principal achado.\n\n"
    "CONTEXTO CIENTÍFICO IMPORTANTE: Este estudo foi criticado por Natália Pasternak "
    "(Instituto Questão de Ciência) e Wellington Zangari (Instituto de Psicologia da USP) por "
    "uma falha fundamental de desenho — usar parentes de primeiro grau como controles infla "
    "falsos positivos, já que eles compartilham 50% do DNA. A definição do fenótipo (mediunidade "
    "autorrelatada) carece de critérios biológicos objetivos. Não existe replicação independente. "
    "Os próprios autores afirmam que replicação é 'necessária e indispensável.' Esta pesquisa "
    "é preliminar e deve ser interpretada como uma exploração inicial intrigante, não como "
    "ciência validada."
)

NEUROIMAGING_CONTEXT = (
    "Grupos de pesquisa brasileiros produziram estudos significativos de neuroimagem de médiuns "
    "espíritas durante estados de transe:\n\n"
    "• Estudos de SPECT (PLoS ONE 2012) de médiuns espíritas brasileiros experientes mostraram "
    "desativação do lobo frontal durante a psicografia (escrita automática) — consistente com "
    "inibição intencional do controle executivo, e não dissociação patológica.\n\n"
    "• Estudos de EEG (SciELO 2016) encontraram maior potência theta e beta durante experiências "
    "sensoriais anômalas vs. controles. Ondas gama e beta distinguiram estados de comunicação "
    "mediúnica de tarefas mentais comuns.\n\n"
    "• Um achado consistente entre os estudos: o transe mediúnico é um estado alterado de "
    "consciência mensurável, distinto da vigília comum — não fingido. Estudos de saúde mental "
    "de médiuns espíritas brasileiros consistentemente encontram taxas abaixo da média de "
    "transtornos psiquiátricos e níveis socioeducacionais acima da média, sugerindo que a "
    "prática mediúnica em contextos religiosos estruturados é psicologicamente adaptativa.\n\n"
    "• O receptor de serotonina 2A (HTR2A) — analisado neste relatório — é o principal "
    "mediador de estados alterados de consciência tanto em contextos farmacológicos (psilocibina) "
    "quanto não farmacológicos (meditação, transe)."
)

CULTURAL_CONTEXT = (
    "O Brasil tem um contexto cultural singularmente favorável para explorar a interseção "
    "entre genética e experiências espirituais:\n\n"
    "• Censo IBGE 2022: 1,84% dos brasileiros se identificam como espíritas (~3,2 milhões), mas "
    "uma população muito maior entre tradições católicas, de Umbanda, Candomblé e evangélicas "
    "se envolve com a mediunidade culturalmente.\n\n"
    "• O Brasil é amplamente descrito como o maior país espírita do mundo, com mais de 12.000 "
    "instituições espíritas.\n\n"
    "• Perfil demográfico espírita: 60,6% feminino, 48% com diploma universitário, 96,6% "
    "com acesso à internet — um segmento de consumidores premium e letrados digitalmente.\n\n"
    "• O centro de pesquisa NUPES (Prof. Alexander Moreira-Almeida, UFJF) trabalha há "
    "décadas para legitimar o estudo científico de experiências espirituais no Brasil.\n\n"
    "• O estudo da USP de 2025 gerou grande cobertura da mídia (CNN Brasil, Estado de Minas), "
    "demonstrando significativo apetite público por conteúdo que une ciência e espiritualidade."
)


# ---------------------------------------------------------------------------
# AVISOS LEGAIS
# ---------------------------------------------------------------------------

DISCLAIMERS = {
    "personality": (
        "Este relatório fornece informações educativas sobre variantes genéticas que pesquisas "
        "científicas associaram a tendências de personalidade. A personalidade é moldada por uma "
        "interação complexa de fatores genéticos, ambientais, culturais e experienciais. "
        "Estes resultados descrevem tendências probabilísticas baseadas em pesquisas no nível "
        "populacional — eles não determinam quem você é ou como você vai se comportar. Os tamanhos "
        "de efeito são geralmente pequenos, e variantes genéticas explicam apenas uma fração da "
        "variação nos traços de personalidade. Isto não é um diagnóstico médico ou psicológico."
    ),
    "mentalHealth": (
        "Este relatório descreve variantes genéticas que pesquisas científicas associaram a "
        "fatores de risco para certas condições relacionadas ao humor e à ansiedade. Ter uma "
        "variante associada a risco elevado NÃO significa que você desenvolverá qualquer "
        "condição. A maioria das pessoas com essas variantes nunca desenvolve a condição "
        "associada. Este relatório NÃO é uma ferramenta de diagnóstico e não substitui a "
        "avaliação por um profissional de saúde mental qualificado. Se você tem preocupações "
        "sobre sua saúde mental, consulte um profissional de saúde. Este relatório é fornecido "
        "apenas para fins educativos e não constitui aconselhamento médico. "
        "CVV - Centro de Valorização da Vida: ligue 188 (24h)."
    ),
    "spiritualSensitivity": (
        "Este relatório explora pesquisas científicas emergentes sobre a genética de experiências "
        "espirituais, absorção e estados transcendentes. O estudo USP/BJPsych de 2025 "
        "(PMID 39874024) identificou genes candidatos potencialmente associados à mediunidade "
        "em praticantes espíritas — esta pesquisa é preliminar e não foi replicada "
        "independentemente. Nenhum teste genético validado para mediunidade ou habilidade espiritual "
        "existe. Este relatório reflete hipóteses científicas atuais e deve ser interpretado "
        "como uma exploração cultural e científica, não como uma ferramenta preditiva ou "
        "diagnóstica. O GeneHealth respeita todas as tradições espirituais e apresenta estas "
        "informações no espírito de curiosidade científica e apreciação cultural."
    ),
}


# ---------------------------------------------------------------------------
# FUNÇÕES DE ANÁLISE
# ---------------------------------------------------------------------------

def _lookup_snp(variants: Dict[str, Tuple[str, str, str]], rsid: str) -> Optional[str]:
    """Look up a SNP genotype from the variants dict.

    Returns the genotype string (e.g., 'AG') or None if not found.
    The variants dict maps rsid -> (chromosome, position, genotype).
    """
    if rsid in variants:
        return variants[rsid][2]  # genotype is the third element
    return None


def _analyze_snp_category(
    variants: Dict[str, Tuple[str, str, str]],
    snp_database: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """Analyze a category of SNPs against user variants."""
    results = []

    for snp_info in snp_database:
        rsid = snp_info["rsid"]
        genotype = _lookup_snp(variants, rsid)

        if genotype is None:
            # SNP not found in user data
            results.append({
                "name": snp_info["name"],
                "gene": snp_info["gene"],
                "rsid": rsid,
                "genotype": "Não disponível",
                "interpretation": (
                    f"Este SNP ({rsid}) não foi encontrado no seu arquivo de DNA. "
                    f"Isso pode ocorrer porque o chip de genotipagem do seu provedor não "
                    f"inclui esta variante."
                ),
                "effectSize": snp_info["effectSize"],
                "scientificBasis": snp_info["scientificBasis"],
                "references": snp_info["references"],
                "actionableInsights": [],
                "score": None,
            })
            continue

        # Normalize genotype for lookup (try both orientations)
        genotype_upper = genotype.upper().strip()
        reversed_gt = genotype_upper[::-1] if len(genotype_upper) == 2 else genotype_upper

        genotype_info = snp_info["genotypes"].get(genotype_upper)
        if genotype_info is None:
            genotype_info = snp_info["genotypes"].get(reversed_gt)

        if genotype_info is None:
            # Unknown genotype combination
            results.append({
                "name": snp_info["name"],
                "gene": snp_info["gene"],
                "rsid": rsid,
                "genotype": genotype_upper,
                "interpretation": (
                    f"Seu genótipo ({genotype_upper}) para {snp_info['gene']} {rsid} "
                    f"não está no nosso banco de interpretações padrão. Isso pode representar "
                    f"uma variante rara ou uma diferença de orientação de fita."
                ),
                "effectSize": snp_info["effectSize"],
                "scientificBasis": snp_info["scientificBasis"],
                "references": snp_info["references"],
                "actionableInsights": [],
                "score": None,
            })
            continue

        results.append({
            "name": snp_info["name"],
            "gene": snp_info["gene"],
            "rsid": rsid,
            "genotype": genotype_upper,
            "genotypeLabel": genotype_info.get("label", genotype_upper),
            "genotypeNickname": genotype_info.get("nickname", ""),
            "interpretation": genotype_info["interpretation"],
            "effectSize": snp_info["effectSize"],
            "scientificBasis": snp_info["scientificBasis"],
            "references": snp_info["references"],
            "actionableInsights": snp_info["actionableInsights"],
            "score": genotype_info.get("score"),
        })

    return results


def _generate_category_summary(traits: List[Dict[str, Any]], category: str) -> str:
    """Generate a summary paragraph for a category based on analyzed traits."""
    analyzed = [t for t in traits if t["genotype"] != "Não disponível" and t.get("score") is not None]

    if not analyzed:
        return f"Nenhum SNP de {category} foi encontrado nos seus dados de DNA. Seu chip de genotipagem pode não incluir essas variantes."

    avg_score = sum(t["score"] for t in analyzed) / len(analyzed)
    count = len(analyzed)
    total = len(traits)

    if category == "personality":
        if avg_score >= 0.6:
            return (
                f"Com base em {count} de {total} variantes de personalidade analisadas, seu perfil "
                f"genético sugere tendências elevadas para sensibilidade emocional e tomada de risco. "
                f"Lembre-se de que a personalidade é moldada principalmente pela experiência de vida, "
                f"cultura e escolha pessoal — a genética fornece apenas uma pequena parte do quadro."
            )
        elif avg_score >= 0.4:
            return (
                f"Com base em {count} de {total} variantes de personalidade analisadas, seu perfil "
                f"genético mostra uma combinação equilibrada de traços. Você carrega uma combinação "
                f"típica de variantes relacionadas à personalidade. A personalidade é predominantemente "
                f"moldada pelo ambiente e experiência."
            )
        else:
            return (
                f"Com base em {count} de {total} variantes de personalidade analisadas, seu perfil "
                f"genético sugere regulação emocional típica e tolerância moderada ao risco. "
                f"Estas são tendências no nível populacional — sua personalidade real é moldada "
                f"principalmente por suas experiências de vida únicas."
            )

    elif category == "mentalHealth":
        if avg_score >= 0.6:
            return (
                f"Com base em {count} de {total} variantes de resposta ao estresse analisadas, seu "
                f"perfil genético mostra algumas variantes associadas a sensibilidade aumentada ao "
                f"estresse. Isso é apenas informativo — a maioria das pessoas com essas variantes "
                f"nunca desenvolve condições de saúde mental. Manejo proativo do estresse, exercício "
                f"regular e forte apoio social são benéficos para todos. Consulte um profissional "
                f"de saúde para qualquer preocupação de saúde mental."
            )
        elif avg_score >= 0.4:
            return (
                f"Com base em {count} de {total} variantes analisadas, seu perfil genético de "
                f"resposta ao estresse está dentro de faixas típicas. Você carrega uma combinação "
                f"comum de variantes protetoras e de sensibilidade. A saúde mental é moldada por "
                f"muitos fatores além da genética."
            )
        else:
            return (
                f"Com base em {count} de {total} variantes analisadas, seu perfil genético "
                f"mostra variantes predominantemente protetoras para resposta ao estresse. Isso é "
                f"encorajador, mas não garante resultados de saúde mental — ambiente, estilo de "
                f"vida e eventos da vida desempenham papéis importantes."
            )

    else:  # spiritualSensitivity
        htr2a_traits = [t for t in analyzed if t["gene"] == "HTR2A"]
        htr2a_score = sum(t["score"] for t in htr2a_traits) / len(htr2a_traits) if htr2a_traits else 0.5

        if htr2a_score >= 0.6:
            return (
                f"Com base em {count} de {total} variantes analisadas, seu perfil serotoninérgico "
                f"sugere sensibilidade potencialmente aumentada a experiências transcendentes e "
                f"contemplativas. Suas variantes de HTR2A (receptor de serotonina 2A) estão associadas "
                f"a densidade alterada de receptores — pesquisas mostram que isso correlaciona com a "
                f"intensidade de experiências místicas. Esta é uma tendência genética, não um diagnóstico "
                f"espiritual. O estudo da USP de 2025 explorou temas semelhantes através de uma lente "
                f"diferente (genética do sistema imune em médiuns), e ambas as linhas de pesquisa "
                f"apontam para a biologia complexa subjacente às experiências espirituais."
            )
        elif htr2a_score >= 0.4:
            return (
                f"Com base em {count} de {total} variantes analisadas, seu perfil serotoninérgico "
                f"mostra um padrão de sensibilidade intermediário. Suas variantes de HTR2A sugerem "
                f"densidade e sensibilidade típicas de receptores a estados alterados. Experiências "
                f"espirituais são moldadas por prática, cultura, crença e biologia em conjunto."
            )
        else:
            return (
                f"Com base em {count} de {total} variantes analisadas, seu perfil serotoninérgico "
                f"mostra um padrão de sensibilidade padrão. Experiências espirituais e contemplativas "
                f"são acessíveis a todas as pessoas independentemente do perfil genético — prática e "
                f"intenção são os principais impulsionadores."
            )


def _generate_overall_profile(traits: List[Dict[str, Any]], category: str) -> str:
    """Generate a one-line overall profile label."""
    analyzed = [t for t in traits if t.get("score") is not None]
    if not analyzed:
        return "Dados insuficientes"

    avg = sum(t["score"] for t in analyzed) / len(analyzed)

    profiles = {
        "personality": {
            0.6: "Sensível e Explorador",
            0.4: "Equilibrado e Adaptável",
            0.0: "Estável e Comedido",
        },
        "mentalHealth": {
            0.6: "Sensibilidade Aumentada ao Estresse",
            0.4: "Resposta ao Estresse Típica",
            0.0: "Perfil Resiliente",
        },
        "spiritualSensitivity": {
            0.6: "Sensibilidade Transcendente Elevada",
            0.4: "Sensibilidade Intermediária",
            0.0: "Sensibilidade Padrão",
        },
    }

    cat_profiles = profiles.get(category, profiles["personality"])
    for threshold in sorted(cat_profiles.keys(), reverse=True):
        if avg >= threshold:
            return cat_profiles[threshold]
    return "Típico"


# ---------------------------------------------------------------------------
# FUNÇÃO PRINCIPAL DE ANÁLISE
# ---------------------------------------------------------------------------

def analyze_mind_spirit(
    variants: Dict[str, Tuple[str, str, str]]
) -> Dict[str, Any]:
    """
    Run full Mind & Spirit analysis on user variants.

    Args:
        variants: Dict mapping rsid -> (chromosome, position, genotype)

    Returns:
        Dict with complete analysis results
    """
    # Analyze each category
    personality_traits = _analyze_snp_category(variants, PERSONALITY_SNPS)
    mental_health_traits = _analyze_snp_category(variants, MENTAL_HEALTH_SNPS)
    spiritual_traits = _analyze_snp_category(variants, SPIRITUAL_SNPS)

    # Also include OXTR rs53576 in spiritual context (absorption/empathy)
    # It's already analyzed in personality — reference it
    oxtr_in_spiritual = None
    for t in personality_traits:
        if t["rsid"] == "rs53576":
            oxtr_in_spiritual = {
                **t,
                "name": "Absorção e Sensibilidade Empática (Contexto Espiritual)",
                "scientificBasis": (
                    t["scientificBasis"] + " No contexto de experiências espirituais, "
                    "variantes de OXTR foram associadas ao constructo psicológico de "
                    "'absorção' — a tendência de se tornar profundamente imerso em experiências. "
                    "A absorção tem ~50% de herdabilidade em estudos com gêmeos e é o traço de "
                    "personalidade mais consistentemente associado a experiências mediúnicas e "
                    "místicas relatadas em diversas culturas."
                ),
            }
            break

    if oxtr_in_spiritual:
        spiritual_traits.append(oxtr_in_spiritual)

    return {
        "personality_traits": personality_traits,
        "mental_health_traits": mental_health_traits,
        "spiritual_traits": spiritual_traits,
    }


def generate_mind_spirit_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate the final JSON output for the Mind & Spirit report.

    Args:
        result: Output from analyze_mind_spirit()

    Returns:
        Dict matching the MindSpiritReport frontend JSON schema
    """
    personality_traits = result["personality_traits"]
    mental_health_traits = result["mental_health_traits"]
    spiritual_traits = result["spiritual_traits"]

    return {
        "reportType": "mind_spirit",
        "version": "1.0",
        "personality": {
            "summary": _generate_category_summary(personality_traits, "personality"),
            "overallProfile": _generate_overall_profile(personality_traits, "personality"),
            "traits": personality_traits,
        },
        "mentalHealth": {
            "summary": _generate_category_summary(mental_health_traits, "mentalHealth"),
            "overallProfile": _generate_overall_profile(mental_health_traits, "mentalHealth"),
            "disclaimer": DISCLAIMERS["mentalHealth"],
            "traits": mental_health_traits,
        },
        "spiritualSensitivity": {
            "summary": _generate_category_summary(spiritual_traits, "spiritualSensitivity"),
            "overallProfile": _generate_overall_profile(spiritual_traits, "spiritualSensitivity"),
            "uspStudyContext": USP_STUDY_CONTEXT,
            "neuroimagingContext": NEUROIMAGING_CONTEXT,
            "culturalContext": CULTURAL_CONTEXT,
            "traits": spiritual_traits,
        },
        "disclaimers": DISCLAIMERS,
    }
