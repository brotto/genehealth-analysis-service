"""
Analisador de Sonhos & Arquitetura do Sono
Analisa variantes genéticas associadas à vivacidade dos sonhos, arquitetura
do sono REM, cronótipo circadiano e fenômenos noturnos (paralisia do sono,
hipnagogia, propensão a sonhos lúcidos).

QUATRO CATEGORIAS:
  1. Vivacidade & Memória dos Sonhos  — HTR2A, DRD4, FAAH, BDNF (cross-ref)
  2. Arquitetura do Sono               — PER3, ADORA2A, NRXN1
  3. Cronótipo & Ritmo Circadiano      — CLOCK, CRY1, PER1, PER2, ARNTL
  4. Fenômenos Noturnos                — HTR2A (contexto paralisia/hipnagogia),
                                         COMT (cross-ref), DRD4 (sonhos lúcidos)

BASE CIENTÍFICA:
  - Sonhos & HTR2A: Nichols (2004) Pharmacol Ther — receptores 5-HT2A e conteúdo onírico
  - HTR2A rs6311: Berger et al. (2003) Neuropsychobiology — densidade de receptores e REM
  - DRD4 rs1800955: Schmack et al. (2013) PLOS ONE — dopamina e conteúdo dos sonhos
  - FAAH rs324420: Fezza et al. (2014) Molecules — endocannabinoides e sono REM
  - PER3 VNTR / rs228697: Viola et al. (2007) Current Biology — arquitetura do sono
  - ADORA2A rs5751876: Rétey et al. (2005) PNAS — pressão do sono e adenosina
  - CLOCK rs1801260: Mishima et al. (2005) Sleep — ritmo circadiano e atraso de fase
  - CRY1 rs2287161: Patke et al. (2017) Cell — período circadiano +~50 min (DSPD)
  - PER1 rs2735611: Katzenberg et al. (1998) / Carpen et al. (2006) — matutinidade
  - PER2 rs2304672: Archer et al. (2003) Sleep — cronótipo vespertino/matutino

CONEXÃO ESPIRITUAL (identidade GeneHealth):
  - HTR2A está no centro tanto dos sonhos vívidos quanto da sensibilidade espiritual
  - Paralisia do sono e hipnagogia são frequentemente reportadas como experiências
    espirituais em diversas culturas (visões, sensação de presença, viagens astrais)
  - Sonhos lúcidos têm base neurobiológica em dopamina prefrontal (DRD4, COMT)

RESSALVAS:
  - Todas as associações são probabilísticas (cada SNP explica <2% da variância)
  - Ambiente, higiene do sono e substâncias modificam fortemente qualquer tendência
  - Não diagnóstico — para distúrbios do sono, consulte um especialista em medicina do sono
"""

from typing import Dict, List, Tuple, Any, Optional


# ---------------------------------------------------------------------------
# UTILITÁRIOS
# ---------------------------------------------------------------------------

COMPLEMENT: Dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _complement(allele: str) -> str:
    return COMPLEMENT.get(allele.upper(), allele.upper())


def _count_allele(genotype: str, allele: str) -> int:
    """Conta ocorrências de um alelo no genótipo com handling de strand flip."""
    genotype = genotype.upper().replace("-", "")
    allele = allele.upper()
    comp = _complement(allele)
    count = genotype.count(allele)
    if count == 0 and comp != allele:
        count = genotype.count(comp)
    return min(count, 2)


def _analyze_snp_list(
    variants: Dict[str, Tuple[str, str, str]],
    snp_list: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    Processa uma lista de SNPs contra o dicionário de variantes do usuário.
    Retorna lista de resultados enriquecidos com genótipo e interpretação.
    """
    results = []
    for snp in snp_list:
        rsid = snp["rsid"].lower()
        chrom, pos, genotype = variants.get(rsid, (None, None, None))

        if genotype:
            genotype = genotype.upper().replace("-", "")
            geno_data = snp["genotypes"].get(genotype)

            # Tenta complemento se não encontrado diretamente
            if not geno_data:
                comp_geno = "".join(_complement(a) for a in genotype)
                geno_data = snp["genotypes"].get(comp_geno)
                # Tenta ordem inversa
                if not geno_data and len(genotype) == 2:
                    geno_data = snp["genotypes"].get(genotype[::-1])
                if not geno_data:
                    rev_comp = comp_geno[::-1]
                    geno_data = snp["genotypes"].get(rev_comp)

            if not geno_data:
                # Fallback: genótipo encontrado mas não mapeado
                geno_data = {
                    "label": genotype,
                    "nickname": "Variante não catalogada",
                    "interpretation": (
                        f"Seu genótipo {genotype} foi detectado mas não está "
                        "catalogado neste SNP. Pode ser um alelo raro."
                    ),
                    "score": 0.3,
                }

            results.append({
                **snp,
                "genotype": genotype,
                "genotypeData": geno_data,
                "found": True,
            })
        else:
            results.append({
                **snp,
                "genotype": "Not available",
                "genotypeData": {
                    "label": "N/D",
                    "nickname": "Não disponível",
                    "interpretation": (
                        "Este SNP não foi encontrado no seu arquivo genômico. "
                        "Kits de diferentes empresas cobrem conjuntos distintos de variantes."
                    ),
                    "score": None,
                },
                "found": False,
            })

        results[-1]["category"] = snp.get("category", "")

    return results


# ---------------------------------------------------------------------------
# BANCO DE SNPs: Vivacidade & Memória dos Sonhos
# ---------------------------------------------------------------------------

DREAM_VIVIDNESS_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs6311",
        "gene": "HTR2A",
        "name": "Vivacidade dos Sonhos e Densidade REM",
        "chromosome": "13",
        "position": 46897343,
        "effectSize": "moderate",
        "category": "dreamVividness",
        "genotypes": {
            "GG": {
                "label": "G/G (−1438G)",
                "nickname": "Sonhador Intenso",
                "interpretation": (
                    "Duas cópias do alelo G associado a maior expressão do receptor 5-HT2A no "
                    "córtex. Receptores 5-HT2A mais densos estão vinculados a sonhos mais vívidos, "
                    "emocionalmente intensos e memoráveis — o mesmo mecanismo que explica os efeitos "
                    "oníricos de substâncias psicodélicas (que ativam exatamente este receptor). "
                    "Você pode ter uma vida onírica naturalmente rica e intensa."
                ),
                "score": 0.8,
            },
            "AG": {
                "label": "A/G",
                "nickname": "Sonhador Moderado",
                "interpretation": (
                    "Genótipo heterozigoto para o promotor do HTR2A. Expressão de 5-HT2A "
                    "intermediária. Você provavelmente experimenta sonhos com intensidade variável "
                    "— às vezes muito vívidos, outras vezes mais difusos. A qualidade dos sonhos "
                    "pode ser sensível a fatores como estresse, alimentação e uso de substâncias."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "A/A (−1438A)",
                "nickname": "Sonhador Sereno",
                "interpretation": (
                    "Duas cópias do alelo A. Expressão levemente reduzida do receptor 5-HT2A. "
                    "Seus sonhos tendem a ser menos intrusivos emocionalmente. Você pode lembrar "
                    "menos dos sonhos ao acordar, mas pode desfrutar de um sono REM mais tranquilo "
                    "e reparador."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "HTR2A rs6311 (C−1438T/A) está no promotor do gene do receptor de serotonina 2A. "
            "O alelo G aumenta a transcrição de HTR2A em ~30% comparado ao alelo A (Polesskaya & "
            "Sokolov, 2002). No sono, os receptores 5-HT2A são ativos durante a vigília e REM, "
            "regulando a intensidade e o conteúdo emocional dos sonhos. Berger et al. (2003) "
            "(Neuropsychobiology 47:2) demonstraram correlação entre variantes promotoras do HTR2A "
            "e a arquitetura do sono REM. A conexão é mecanisticamente clara: psicodélicos clássicos "
            "(LSD, psilocibina) que produzem experiências oníricas intensas agem como agonistas "
            "exatamente neste receptor. Este SNP também está presente no mind_spirit_analyzer da "
            "GeneHealth no contexto de sensibilidade espiritual — a mesma variante que amplia a "
            "experiência onírica tende a ampliar a sensibilidade transcendente na vigília."
        ),
        "references": [
            {"pmid": "12566938", "title": "HTR2A promoter polymorphism and receptor expression", "year": 2003},
            {"pmid": "15163253", "title": "Serotonin 2A receptor and sleep architecture", "year": 2004},
        ],
        "actionableInsights": [
            "Portadores G/G: Evite telas brilhantes 2h antes de dormir — a estimulação serotonérgica pode intensificar demais os sonhos",
            "Meditação noturna pode ajudar a 'preparar' o espaço mental para sonhos mais integrativos",
            "Se os sonhos forem perturbadores, estratégias de processamento onírico (journaling) são especialmente valiosas para você",
            "Portadores A/A: Técnicas de indução de sonhos lúcidos (MILD, WILD) podem ser necessárias para aumentar a recall de sonhos",
        ],
    },
    {
        "rsid": "rs1800955",
        "gene": "DRD4",
        "name": "Dopamina, Novidade Onírica e Sonhos Lúcidos",
        "chromosome": "11",
        "position": 637304,
        "effectSize": "moderate",
        "category": "dreamVividness",
        "genotypes": {
            "TT": {
                "label": "T/T",
                "nickname": "Conteúdo Onírico Típico",
                "interpretation": (
                    "Genótipo de referência para DRD4 −521. Expressão e função típica do receptor "
                    "de dopamina D4. Seus sonhos tendem a refletir predominantemente conteúdo do "
                    "cotidiano e memórias recentes."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Sonhador Explorador",
                "interpretation": (
                    "Uma cópia do alelo C. O receptor D4 com variante C tem expressão reduzida no "
                    "córtex pré-frontal, o que pode aumentar a probabilidade de conteúdo onírico "
                    "incomum, bizarro ou de exploração de novidades. Alguns estudos associam "
                    "variantes DRD4 a propensão para sonhos lúcidos."
                ),
                "score": 0.6,
            },
            "CC": {
                "label": "C/C",
                "nickname": "Sonhador Visionário",
                "interpretation": (
                    "Duas cópias do alelo C. Variante associada ao receptor D4 de menor afinidade "
                    "pela dopamina. Correlacionada com busca por novidades, criatividade e, em "
                    "contexto onírico, conteúdo mais incomum e potencialmente mais bizarro ou "
                    "visionário. A atividade dopaminérgica prefrontal durante o REM está ligada "
                    "à metacognição onírica (consciência de que se está sonhando)."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "DRD4 −521C/T (rs1800955) afeta a expressão do receptor de dopamina D4 — variantes de "
            "menor expressão no córtex pré-frontal foram associadas à busca por novidades (Strobel "
            "et al., 1999). Durante o sono REM, a dopamina prefrontal é crucial para a metacognição "
            "onírica. Schmack et al. (2013, PLOS ONE) demonstraram que perfis dopaminérgicos "
            "influenciam o conteúdo e a bizarria dos sonhos. A ligação com sonhos lúcidos é explorada "
            "por Hobson & Friston (2012): a consciência durante o REM requer ativação de áreas "
            "frontais que normalmente estão inativas — padrão favorecido por certos perfis de "
            "receptor D4. O DRD4 de 7 repetições (VNTR, bem correlacionado com rs1800955) é o alelo "
            "mais fortemente associado à novidade, criatividade e comportamentos exploratórios."
        ),
        "references": [
            {"pmid": "23516459", "title": "Dopaminergic prediction errors drive dream content", "year": 2013},
            {"pmid": "10205490", "title": "DRD4 promoter polymorphism and novelty seeking", "year": 1999},
        ],
        "actionableInsights": [
            "Portadores CC: Técnicas de realidade verificada (reality checks) durante o dia podem aumentar a probabilidade de sonhos lúcidos",
            "Mantenha um diário de sonhos ao lado da cama — sua predisposição a conteúdo incomum vale registrar",
            "Portadores TT: Técnicas de MILD (Mnemonic Induction of Lucid Dreams) podem ser especialmente eficazes por fortalecer a intenção pré-sono",
        ],
    },
    {
        "rsid": "rs324420",
        "gene": "FAAH",
        "name": "Endocannabinoides, Intensidade Emocional dos Sonhos",
        "chromosome": "1",
        "position": 46891060,
        "effectSize": "moderate",
        "category": "dreamVividness",
        "genotypes": {
            "AA": {
                "label": "C385A / A385A",
                "nickname": "Sonhos Emocionalmente Amplificados",
                "interpretation": (
                    "Portador do alelo A (385A). A enzima FAAH com 385A tem atividade reduzida, "
                    "resultando em níveis mais altos de anandamida (o 'canabinóide endógeno' ou "
                    "'molécula da bem-aventurança'). Maior anandamida está associada a sonhos mais "
                    "intensos emocionalmente, com maior tonalidade positiva e, em alguns estudos, "
                    "menor ansiedade onírica. Você pode ter uma vida de sonhos emocionalmente "
                    "rica com tonalidade hedônica mais elevada."
                ),
                "score": 0.7,
            },
            "AC": {
                "label": "Heterozigoto C385A",
                "nickname": "Modulação Endocannabinoide Balanceada",
                "interpretation": (
                    "Uma cópia de cada variante. Atividade intermediária da FAAH e níveis de "
                    "anandamida entre os extremos. Você provavelmente experimenta uma vida onírica "
                    "com intensidade emocional variável conforme seu estado emocional da vigília."
                ),
                "score": 0.5,
            },
            "CC": {
                "label": "C/C (385Pro)",
                "nickname": "Processamento Endocannabinoide Padrão",
                "interpretation": (
                    "Genótipo de referência com atividade normal da FAAH. Níveis padrão de "
                    "anandamida. Seus sonhos tendem a ter intensidade emocional típica, sem "
                    "amplificação endocannabinoide particular."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "FAAH (Fatty Acid Amide Hydrolase) rs324420 (C385A) é a variante mais estudada no "
            "sistema endocannabinoide. O alelo 385A produz uma proteína ~30% menos estável, "
            "resultando em menor degradação de anandamida. Prather et al. (2013, PNAS) mostraram "
            "que portadores do alelo A têm resposta reduzida ao estresse e humor mais positivo. "
            "No contexto do sono, o sistema endocannabinoide regula o sono REM — animais sem FAAH "
            "(e portanto com alta anandamida) mostram aumento de sono REM com maior atividade "
            "neuronal durante essa fase. Bhattacharya et al. (2021, Int J Mol Sci) revisaram a "
            "relação entre endocannabinoides e regulação do sono, confirmando o papel da anandamida "
            "na intensidade do sono REM. Este SNP é particularmente relevante para o perfil "
            "espiritual da GeneHealth: a anandamida é chamada de 'molécula da bliss' e sua "
            "conexão com experiências de absorção e bem-estar espiritual está documentada."
        ),
        "references": [
            {"pmid": "23754379", "title": "FAAH C385A and stress reactivity in humans", "year": 2013},
            {"pmid": "34444209", "title": "Endocannabinoid system and sleep regulation", "year": 2021},
        ],
        "actionableInsights": [
            "Portadores A/A: Práticas de mindfulness que aumentam anandamida endógena (meditação, exercício) podem potencializar sua experiência onírica",
            "Portadores C/C: O exercício físico eleva temporariamente a anandamida e pode enriquecer seus sonhos nas noites seguintes",
            "Todos os genótipos: Dieta rica em ácidos graxos ômega-3 suporta a síntese de endocannabinoides",
        ],
    },
]


# ---------------------------------------------------------------------------
# BANCO DE SNPs: Arquitetura do Sono
# ---------------------------------------------------------------------------

SLEEP_ARCHITECTURE_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs228697",
        "gene": "PER3",
        "name": "Arquitetura do Sono REM e Pressão do Sono",
        "chromosome": "1",
        "position": 7847208,
        "effectSize": "strong",
        "category": "sleepArchitecture",
        "genotypes": {
            "CC": {
                "label": "Alelo curto (proxy 4-repeat)",
                "nickname": "Sono Eficiente",
                "interpretation": (
                    "Associado ao alelo PER3 de 4 repetições (proxy). Portadores tendem a acumular "
                    "pressão de sono mais rapidamente, adormecer com facilidade e ter sono profundo "
                    "eficiente. No entanto, são mais sensíveis aos efeitos cognitivos da privação de "
                    "sono. O sono REM é tipicamente mais compacto e concentrado no final da noite."
                ),
                "score": 0.4,
            },
            "CT": {
                "label": "Heterozigoto",
                "nickname": "Perfil de Sono Misto",
                "interpretation": (
                    "Genótipo heterozigoto. Combina características de ambos os perfis de "
                    "arquitetura do sono. Você provavelmente experimenta boa eficiência do sono "
                    "com uma vida de sonhos moderadamente rica — o equilíbrio mais comum."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "Alelo longo (proxy 5-repeat)",
                "nickname": "Sonhador REM Rico",
                "interpretation": (
                    "Associado ao alelo PER3 de 5 repetições (proxy). Este perfil é fascinante: "
                    "portadores têm sono REM mais distribuído ao longo da noite, com maior "
                    "densidade de sonhos e melhor recall onírico. São mais resistentes à privação "
                    "de sono cognitivamente (precisam dormir menos para funcionar bem), mas têm "
                    "maior latência de início do sono. Viola et al. (2007) mostraram que portadores "
                    "5/5 têm ondas delta de maior amplitude — sono profundo mais intenso."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "PER3 (Period Circadian Regulator 3) contém um VNTR (rs57875989) de 4 ou 5 repetições "
            "de 54 pares de bases que é o SNP de sono mais estudado na genômica do consumidor. "
            "rs228697 é um proxy de alta correlação disponível em arrays de genotyping padrão. "
            "Viola et al. (2007, Current Biology 17:4) demonstraram que portadores 5/5 têm mais "
            "ondas lentas (SWA) e maior atividade EEG theta durante a vigília após privação do sono. "
            "Dijk & Archer (2010, J Sleep Res) revisaram extensivamente como o genótipo PER3 "
            "modula tanto a homeostase do sono quanto a resposta à privação. A relevância para "
            "a GeneHealth: portadores do alelo longo têm experiências oníricas naturalmente mais "
            "ricas — um dado que se conecta diretamente ao perfil de sensibilidade espiritual."
        ),
        "references": [
            {"pmid": "17320392", "title": "PER3 polymorphism and sleep homeostasis", "year": 2007},
            {"pmid": "20082661", "title": "PER3 VNTR and sleep regulation", "year": 2010},
        ],
        "actionableInsights": [
            "Portadores TT (proxy 5/5): Garanta 8+ horas — sua arquitetura de sono é naturalmente mais longa para completar ciclos REM",
            "Portadores CC (proxy 4/4): Você pode ser mais resistente cognitivamente à privação, mas não negligencie o sono de qualidade",
            "Todos: Manter horários regulares de sono maximiza a densidade do REM independentemente do genótipo",
        ],
    },
    {
        "rsid": "rs5751876",
        "gene": "ADORA2A",
        "name": "Pressão do Sono e Sensibilidade à Cafeína",
        "chromosome": "22",
        "position": 24822676,
        "effectSize": "strong",
        "category": "sleepArchitecture",
        "genotypes": {
            "TT": {
                "label": "T/T (1083T)",
                "nickname": "Alta Pressão de Sono",
                "interpretation": (
                    "Genótipo associado a maior sensibilidade ao receptor de adenosina A2A. "
                    "A adenosina é a principal molécula que acumula a 'pressão de sono' durante "
                    "a vigília. Portadores TT sentem mais sono ao longo do dia, têm maior "
                    "propensão a adormecer e são mais sensíveis à cafeína (que bloqueia exatamente "
                    "este receptor). Quando dormem, tendem a ter sono mais profundo e reparador. "
                    "Seus sonhos podem ser mais intensos após privação parcial."
                ),
                "score": 0.7,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Sensibilidade Intermediária",
                "interpretation": (
                    "Perfil heterozigoto. Sensibilidade à adenosina e à cafeína entre os extremos. "
                    "Você provavelmente tolera cafeína moderada sem grande impacto no sono se "
                    "consumida antes das 14h."
                ),
                "score": 0.5,
            },
            "CC": {
                "label": "C/C (1083C)",
                "nickname": "Sono Resiliente",
                "interpretation": (
                    "Menor sensibilidade ao receptor A2A de adenosina. Portadores CC tendem a "
                    "acumular pressão de sono mais lentamente, são menos sonolentos durante o dia "
                    "e mais tolerantes à cafeína noturna. Podem funcionar melhor com privação "
                    "parcial de sono, mas o sono profundo (ondas lentas) pode ser levemente "
                    "reduzido em comparação com portadores TT."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "ADORA2A rs5751876 é o SNP de sensibilidade à cafeína mais bem validado. Rétey et al. "
            "(2005, PNAS 102:17) mostraram que portadores T/T têm maior sensibilidade subjetiva e "
            "objetiva à cafeína (medida por EEG e escalas de sonolência). O receptor A2A de adenosina "
            "é o principal alvo da cafeína — quando a adenosina se liga, sinaliza cansaço; quando "
            "a cafeína bloqueia, você se sente acordado. Portadores T/T têm um receptor que 'amplifica' "
            "tanto o sinal de sono quanto o efeito de vigília da cafeína. Em termos de sonhos: a "
            "arquitetura do sono dos portadores T/T favorece ciclos REM mais intensos quando o sono "
            "finalmente ocorre, potencialmente amplificando a experiência onírica."
        ),
        "references": [
            {"pmid": "16188225", "title": "ADORA2A polymorphism and caffeine sensitivity", "year": 2005},
            {"pmid": "17592536", "title": "Adenosine A2A receptor and sleep homeostasis", "year": 2007},
        ],
        "actionableInsights": [
            "Portadores T/T: Evite cafeína após as 13h — seu sistema adenosinérgico amplifica o impacto",
            "Portadores C/C: Você pode tolerar café vespertino, mas monitore a qualidade do sono individualmente",
            "Todos: A 'dívida de sono' acumulada leva a sonhos mais intensos na recuperação — use isso conscientemente",
        ],
    },
]


# ---------------------------------------------------------------------------
# BANCO DE SNPs: Cronótipo & Ritmo Circadiano
# ---------------------------------------------------------------------------

CHRONOTYPE_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs1801260",
        "gene": "CLOCK",
        "name": "Ritmo Circadiano e Atraso de Fase do Sono",
        "chromosome": "4",
        "position": 56292869,
        "effectSize": "moderate",
        "category": "chronotype",
        "genotypes": {
            "CC": {
                "label": "C/C (3111C)",
                "nickname": "Ciclo Circadiano Padrão",
                "interpretation": (
                    "Genótipo de referência para CLOCK 3111T/C. Período circadiano típico (~24h). "
                    "Você provavelmente se adapta bem a horários convencionais de sono e tem "
                    "relativa facilidade para ajustar o ciclo sono-vigília."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Tendência Noturna Leve",
                "interpretation": (
                    "Uma cópia do alelo T. Associado a uma leve tendência para cronótipo vespertino "
                    "(noturno). Você pode ter facilidade para ficar acordado até mais tarde e "
                    "dificuldade para acordar cedo, especialmente na adolescência e início da "
                    "vida adulta."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "T/T (3111T)",
                "nickname": "Cronótipo Noturno",
                "interpretation": (
                    "Duas cópias do alelo T. Mishima et al. (2005) associaram este genótipo a "
                    "atraso de fase do sono e padrão circadiano noturno. Portadores TT têm "
                    "maior tendência a dormir tarde, acordar tarde e funcionar melhor à noite. "
                    "O pico de melatonina ocorre mais tarde em comparação com portadores CC. "
                    "Sonhos mais intensos tendem a ocorrer próximos ao horário de despertar "
                    "natural tardio deste cronótipo."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "CLOCK 3111T/C (rs1801260) está na região 3' UTR do gene CLOCK, um dos principais "
            "reguladores do relógio circadiano molecular. O CLOCK forma o heterodímero CLOCK:BMAL1 "
            "que ativa a transcrição dos genes Per e Cry. Mishima et al. (2005, Sleep 28:1) "
            "demonstraram associação significativa do alelo T com síndrome do atraso de fase do "
            "sono (DSPS) em japoneses. Katzenberg et al. (1998, Sleep 21:6) também reportaram "
            "associação com matutinidade vs. vespertinidade. O CLOCK é o 'maestro' do ritmo "
            "circadiano — variantes aqui afetam não apenas o sono, mas também o timing de picos "
            "hormonais, temperatura corporal e, por consequência, quando os sonhos mais vívidos "
            "ocorrem durante a arquitetura de sono personalizada."
        ),
        "references": [
            {"pmid": "15700720", "title": "CLOCK mutation and delayed sleep phase syndrome", "year": 2005},
            {"pmid": "9756555", "title": "CLOCK gene polymorphism and human diurnal preference", "year": 1998},
        ],
        "actionableInsights": [
            "Portadores T/T: Respeite sua biologia — forçar um ciclo muito matutino pode cronicamente reduzir a qualidade do sono REM",
            "Exposição à luz solar matinal ajuda a 'ancorar' o ciclo circadiano independentemente do genótipo",
            "Portadores T/T que precisam de horários matutinos: melatonina de baixa dose (0,5mg) 5h antes do sono desejado pode ajudar",
        ],
    },
    {
        "rsid": "rs2287161",
        "gene": "CRY1",
        "name": "Período Circadiano Longo e Sono Atrasado",
        "chromosome": "12",
        "position": 107484498,
        "effectSize": "strong",
        "category": "chronotype",
        "genotypes": {
            "CC": {
                "label": "C/C",
                "nickname": "Período Circadiano Típico",
                "interpretation": (
                    "Genótipo de referência. Período circadiano endógeno próximo a 24h. "
                    "Você se adapta bem a horários convencionais e tende a ter boa sincronização "
                    "entre seu relógio interno e o ciclo claro-escuro ambiental."
                ),
                "score": 0.3,
            },
            "CG": {
                "label": "C/G",
                "nickname": "Período Levemente Alongado",
                "interpretation": (
                    "Uma cópia do alelo G. Seu período circadiano endógeno pode ser ligeiramente "
                    "mais longo que 24h. Isso cria uma tendência natural a atrasar o ciclo "
                    "sono-vigília progressivamente se não houver âncoras de luz externas fortes."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "G/G",
                "nickname": "Período Circadiano Alongado (~50 min+)",
                "interpretation": (
                    "Variante de alto impacto descrita por Patke et al. (2017, Cell). Portadores "
                    "têm um período circadiano endógeno mais longo — aproximadamente 50 minutos "
                    "acima da média. Isso significa que seu relógio interno prefere 'correr mais "
                    "lento', resultando em tendência forte a dormis tarde e acordar tarde. Em "
                    "ambientes livres (sem obrigações de horário), portadores naturalmente dormem "
                    "2-3h mais tarde que a média. Esta é uma das variantes de cronótipo mais "
                    "potentes identificadas até hoje."
                ),
                "score": 0.85,
            },
        },
        "scientificBasis": (
            "CRY1 (Cryptochrome 1) é um repressor central do relógio circadiano molecular — "
            "ele inibe o heterodímero CLOCK:BMAL1, fechando o loop de feedback. Patke et al. "
            "(2017, Cell 169:2) identificaram uma mutação de splicing em CRY1 associada à "
            "síndrome do atraso de fase do sono familiar (FASP). rs2287161 é um marcador de "
            "desequilíbrio de ligação com esta variante funcional. Portadores com período mais "
            "longo têm mais dificuldade de sincronização com o ciclo solar padrão — sua 'meia-noite "
            "biológica' ocorre mais tarde. Estudos de sleep EEG mostram que portadores com período "
            "CRY1 longo têm maior intensidade de ondas lentas acumulada quando finalmente dormem, "
            "potencialmente resultando em REM de recuperação mais rico."
        ),
        "references": [
            {"pmid": "28388406", "title": "CRY1 splice variant and delayed sleep phase disorder", "year": 2017},
            {"pmid": "30773234", "title": "CRY1 functional variants and circadian period", "year": 2019},
        ],
        "actionableInsights": [
            "Portadores G/G: Seu cronótipo tem uma base genética forte — validar isso pode aliviar a culpa social de ser 'noturno'",
            "Exposição matinal intensa à luz (10.000 lux por 30 min) é a intervenção mais eficaz para avançar o ciclo",
            "Se trabalho permite, considere ajustar horários para respeitar seu pico cognitivo natural (tipicamente 11h-13h para noturnos)",
        ],
    },
    {
        "rsid": "rs2735611",
        "gene": "PER1",
        "name": "Preferência Matutina e Antecipação Circadiana",
        "chromosome": "17",
        "position": 8123949,
        "effectSize": "moderate",
        "category": "chronotype",
        "genotypes": {
            "AA": {
                "label": "A/A",
                "nickname": "Cronótipo Matutino",
                "interpretation": (
                    "Associado a preferência matutina. Portadores tendem a adormecer mais cedo, "
                    "acordar naturalmente cedo e ter pico de energia e cognição nas primeiras horas "
                    "da manhã. O sono REM mais intenso ocorre nos ciclos do início da madrugada, "
                    "o que pode resultar em sonhos mais vívidos nas primeiras horas após dormir."
                ),
                "score": 0.7,
            },
            "AG": {
                "label": "A/G",
                "nickname": "Cronótipo Intermediário",
                "interpretation": (
                    "Genótipo heterozigoto. Tendência ao cronótipo intermediário — você "
                    "provavelmente se adapta razoavelmente bem a diferentes horários, "
                    "com leve preferência por acordar em horários moderados."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "G/G",
                "nickname": "Tendência Vespertina",
                "interpretation": (
                    "Associado a menor preferência matutina. Você pode ter mais facilidade "
                    "para ficar acordado à noite e mais dificuldade para acordar cedo. "
                    "Em combinação com outros marcadores circadianos, pode indicar cronótipo "
                    "vespertino moderado."
                ),
                "score": 0.4,
            },
        },
        "scientificBasis": (
            "PER1 (Period 1) é um dos três genes 'period' do relógio circadiano. rs2735611 "
            "foi associado a cronotipo e preferência matutina em estudos GWAS de grande escala. "
            "Jones et al. (2019, Nature Communications) identificaram 351 loci para cronótipo em "
            "~700.000 indivíduos do UK Biobank, com PER1 entre os mais significativos. Carpen et al. "
            "(2006, J Sleep Res) também reportaram associação de variantes PER1 com matutinidade. "
            "PER1 regula a fase do relógio circadiano junto com PER2 e PER3 — variantes neste gene "
            "afetam a velocidade com que o relógio 'avança' em resposta à luz matinal."
        ),
        "references": [
            {"pmid": "30696823", "title": "GWAS of chronotype in UK Biobank identifies 351 loci", "year": 2019},
            {"pmid": "16685255", "title": "PER1 polymorphism and diurnal preference", "year": 2006},
        ],
        "actionableInsights": [
            "Portadores A/A: Aproveite seu pico matinal — agende trabalho criativo e decisões importantes antes do meio-dia",
            "Portadores G/G: Considere blocos de trabalho profundo no final da tarde/início da noite quando sua cognição estiver no pico",
            "Todos: Consistência no horário de despertar (mesmo nos fins de semana) é o fator mais importante para a qualidade do sono",
        ],
    },
]


# ---------------------------------------------------------------------------
# BANCO DE SNPs: Fenômenos Noturnos
# (Paralisia do sono, hipnagogia, sonhos lúcidos, experiências de fronteira)
# ---------------------------------------------------------------------------

NOCTURNAL_PHENOMENA_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs6313",
        "gene": "HTR2A",
        "name": "Paralisia do Sono e Experiências Hipnagógicas",
        "chromosome": "13",
        "position": 46897065,
        "effectSize": "moderate",
        "category": "nocturnal",
        "genotypes": {
            "TT": {
                "label": "T/T",
                "nickname": "Fronteira Sono-Vigília Definida",
                "interpretation": (
                    "Genótipo de referência. Transição sono-vigília típica. Você provavelmente "
                    "experimenta a fronteira entre o sono e a vigília de forma relativamente clara, "
                    "sem fenômenos de fronteira frequentes."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Fronteira Permeável",
                "interpretation": (
                    "Uma cópia do alelo C. Expressão aumentada de HTR2A pode tornar a fronteira "
                    "sono-vigília mais permeável. Você pode experimentar com maior frequência "
                    "alucinações hipnagógicas (ao adormecer) ou hipnopômpicas (ao acordar), "
                    "paralisia do sono ocasional ou sensações vívidas na transição entre estados. "
                    "Estas experiências, embora às vezes assustadoras, têm base neurobiológica "
                    "clara e são reportadas em diversas tradições como experiências espirituais."
                ),
                "score": 0.6,
            },
            "CC": {
                "label": "C/C",
                "nickname": "Fronteira Muito Permeável",
                "interpretation": (
                    "Maior expressão do receptor 5-HT2A. A fronteira entre REM e vigília pode "
                    "ser especialmente permeável para você. Episódios de paralisia do sono com "
                    "conteúdo alucinatório rico, sonhos lúcidos espontâneos e experiências "
                    "hipnagógicas intensas são mais prováveis com este perfil. Em muitas culturas, "
                    "estas experiências de fronteira sono-vigília são interpretadas como contato "
                    "com outras dimensões, visões ou mensagens espirituais — e você pode ter "
                    "predisposição genética para essas experiências de fronteira."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "HTR2A rs6313 é outro polimorfismo no promotor do receptor de serotonina 2A. "
            "A paralisia do sono (SP) ocorre quando há intrusão do estado REM na vigília — "
            "os mecanismos que paralisam a musculatura durante o REM persistem ao acordar, "
            "frequentemente acompanhados por alucinações hipnopômpicas. O sistema serotonérgico "
            "via HTR2A regula a transição REM-vigília. Jalal & Hinton (2013, Psychol Sci) "
            "documentaram extensivamente as alucinações da SP em 30+ culturas como 'visitantes "
            "noturnos', 'demônios', 'entidades' — a mesma neurobiologia expressa em narrativas "
            "espirituais distintas. Variantes que aumentam HTR2A podem aumentar tanto a "
            "probabilidade quanto a intensidade dessas experiências de fronteira."
        ),
        "references": [
            {"pmid": "23983260", "title": "Sleep paralysis hallucinations across cultures", "year": 2013},
            {"pmid": "11774094", "title": "HTR2A polymorphisms and serotonin receptor expression", "year": 2001},
        ],
        "actionableInsights": [
            "Portadores CC: Se a paralisia do sono for perturbadora, dormir de lado (não de costas) reduz significativamente a frequência",
            "Você pode explorar técnicas de indução de sonhos lúcidos — sua fronteira permeável pode ser uma porta de entrada",
            "Compreender a neurobiologia desses estados pode transformar experiências assustadoras em exploração consciente",
            "Práticas meditativas que trabalham com estados hipnagógicos (yoga nidra, NSDR) são especialmente relevantes para seu perfil",
        ],
    },
    {
        "rsid": "rs6265",
        "gene": "BDNF",
        "name": "Consolidação de Memória Onírica e Neuroplasticidade",
        "chromosome": "11",
        "position": 27658369,
        "effectSize": "moderate",
        "category": "nocturnal",
        "genotypes": {
            "CC": {
                "label": "Val/Val",
                "nickname": "Consolidação Onírica Padrão",
                "interpretation": (
                    "Genótipo de referência. Secreção normal de BDNF dependente de atividade. "
                    "A consolidação de memória durante o sono segue o padrão típico — sonhos "
                    "relacionados a eventos recentes e processamento emocional usual."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "Val/Met",
                "nickname": "Consolidação Onírica Alterada",
                "interpretation": (
                    "Uma cópia do alelo Met. A secreção reduzida de BDNF dependente de atividade "
                    "pode afetar como as memórias emocionais são processadas e consolidadas durante "
                    "o sono. Você pode ter sonhos com replay de eventos emocionais intensos com "
                    "maior frequência — seu cérebro pode precisar de mais 'iterações' de REM para "
                    "processar experiências emocionalmente carregadas."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "Met/Met",
                "nickname": "Processamento Onírico Intensificado",
                "interpretation": (
                    "Duas cópias do alelo Met (Val66Met). A redução na liberação de BDNF "
                    "dependente de atividade afeta a plasticidade hipocampal e a consolidação "
                    "de memória. No contexto onírico, isso pode se manifestar como sonhos mais "
                    "frequentes com processamento de memórias emocionais, sonhos recorrentes "
                    "em períodos de estresse (o cérebro tentando 'resolver' via REM), e "
                    "possivelmente maior ruminação onírica. Exercício regular eleva o BDNF "
                    "e pode melhorar a qualidade do processamento durante o sono."
                ),
                "score": 0.65,
            },
        },
        "scientificBasis": (
            "BDNF Val66Met (rs6265) já está no mind_spirit_analyzer da GeneHealth no contexto "
            "de neuroplasticidade diurna. Aqui, o enfoque é o papel do BDNF na consolidação de "
            "memória durante o sono. Walker & Stickgold (2004, Neuron 44:1) estabeleceram que "
            "o sono REM é crítico para a consolidação de memórias procedurais e emocionais. "
            "BDNF é o principal mediador da plasticidade sináptica durante essa consolidação. "
            "Bhattacharya et al. (2015, Neurobiology of Sleep and Circadian Rhythms) revisaram "
            "como variantes de BDNF afetam a consolidação de memória dependente do sono. "
            "Este SNP aparece em dois contextos na plataforma GeneHealth: durante a vigília "
            "(plasticidade e humor) e durante o sono (consolidação onírica e processamento "
            "emocional noturno) — ilustrando a continuidade mente-sono."
        ),
        "references": [
            {"pmid": "15450160", "title": "Sleep and the price of plasticity", "year": 2004},
            {"pmid": "26779544", "title": "BDNF and sleep-dependent memory consolidation", "year": 2015},
        ],
        "actionableInsights": [
            "Portadores TT: O exercício aeróbico (especialmente manhã/tarde) eleva BDNF e pode melhorar dramaticamente a qualidade do processamento onírico",
            "Journaling pré-sono ajuda a 'pré-digerir' conteúdo emocional antes que o REM precise processá-lo",
            "Portadores Val/Val: Seu processamento de memória durante o sono é eficiente — aproveite dormindo bem após aprendizado importante",
        ],
    },
]


# ---------------------------------------------------------------------------
# CONSTANTES DE CONTEXTO (para o relatório frontend)
# ---------------------------------------------------------------------------

SPIRITUAL_DREAM_CONTEXT = (
    "Na GeneHealth, entendemos que a experiência onírica é uma janela única para a "
    "consciência. Culturas de todo o mundo trataram os sonhos como portais para o sagrado "
    "— de incubação onírica nos templos de Asclépio na Grécia antiga, às tradições "
    "xamânicas siberianas, aos sonhos proféticos do Candomblé e das tradições indígenas "
    "brasileiras. A neurociência moderna está começando a compreender por que: o estado "
    "REM compartilha características de estados meditativos profundos e experiências com "
    "substâncias psicodélicas, incluindo ativação das mesmas redes neurais e redução do "
    "Default Mode Network. Seus marcadores genéticos neste relatório refletem sua "
    "predisposição biológica para essa dimensão da experiência humana."
)

DISCLAIMERS = {
    "general": (
        "Este relatório é educativo e não constitui diagnóstico médico. "
        "Genética é apenas um fator — higiene do sono, estresse, medicamentos "
        "e ambiente têm impacto igual ou maior na qualidade do sono e dos sonhos. "
        "Para distúrbios do sono, consulte um especialista em medicina do sono."
    ),
    "chronotype": (
        "Cronótipo tem forte base genética (~50% herdabilidade) mas é moldável "
        "por exposição à luz, horários sociais e rotina. Não use estes resultados "
        "para justificar hábitos prejudiciais à saúde."
    ),
    "nocturnal": (
        "Paralisia do sono e experiências hipnagógicas são fenômenos neurológicos benignos, "
        "presentes em 20-40% da população. Se frequentes ou perturbadores, consulte um "
        "especialista — tratamentos comportamentais simples são altamente eficazes."
    ),
}


# ---------------------------------------------------------------------------
# FUNÇÕES AUXILIARES DE RELATÓRIO
# ---------------------------------------------------------------------------

def _generate_category_summary(traits: List[Dict[str, Any]], category: str) -> str:
    """Gera um sumário textual para uma categoria de traços."""
    found = [t for t in traits if t.get("found", False)]
    if not found:
        return "Nenhuma variante desta categoria encontrada no seu arquivo genômico."

    scores = [
        t["genotypeData"]["score"]
        for t in found
        if t["genotypeData"].get("score") is not None
    ]
    if not scores:
        return "Variantes encontradas, mas sem dados de interpretação suficientes."

    avg = sum(scores) / len(scores)

    summaries = {
        "dreamVividness": {
            0.65: "Seu perfil genético sugere uma vida onírica naturalmente intensa e rica — sonhos vívidos, emocionalmente carregados e com boa memorabilidade são características prováveis.",
            0.45: "Você tem um perfil onírico equilibrado — com potencial para sonhos vívidos em condições favoráveis, especialmente após dias emocionalmente intensos.",
            0.0: "Seu perfil sugere uma abordagem mais serena à vida onírica — sonhos tendem a ser menos intrusivos, com recall moderado.",
        },
        "sleepArchitecture": {
            0.65: "Sua arquitetura de sono tende para ciclos REM ricos e sono profundo intenso — você provavelmente precisa de horas adequadas para completar seus ciclos naturalmente mais longos.",
            0.45: "Arquitetura de sono equilibrada com boa resposta à regularidade de horários.",
            0.0: "Sono eficiente com acúmulo rápido de pressão — você pode funcionar bem com ciclos mais compactos.",
        },
        "chronotype": {
            0.65: "Seu relógio biológico tem forte tendência vespertina (noturna) com base genética sólida.",
            0.45: "Cronótipo intermediário — você se adapta a diferentes horários com relativa facilidade.",
            0.0: "Forte predisposição ao cronótipo matutino — você provavelmente funciona melhor cedo.",
        },
        "nocturnal": {
            0.65: "Alta probabilidade de experiências de fronteira sono-vigília (hipnagogia, paralisia do sono, sonhos lúcidos espontâneos).",
            0.45: "Propensão moderada a fenômenos noturnos — experiências de fronteira são possíveis em condições de privação ou estresse.",
            0.0: "Perfil com transições sono-vigília tipicamente bem definidas.",
        },
    }

    cat_summaries = summaries.get(category, summaries["dreamVividness"])
    for threshold in sorted(cat_summaries.keys(), reverse=True):
        if avg >= threshold:
            return cat_summaries[threshold]
    return "Perfil típico."


def _generate_overall_profile(traits: List[Dict[str, Any]], category: str) -> str:
    """Gera o rótulo de perfil geral."""
    found = [t for t in traits if t.get("found", False)]
    if not found:
        return "Dados insuficientes"

    scores = [
        t["genotypeData"]["score"]
        for t in found
        if t["genotypeData"].get("score") is not None
    ]
    if not scores:
        return "Perfil típico"

    avg = sum(scores) / len(scores)

    profiles = {
        "dreamVividness": {
            0.65: "Sonhador Intenso",
            0.45: "Sonhador Equilibrado",
            0.0: "Sonhador Sereno",
        },
        "sleepArchitecture": {
            0.65: "Arquitetura REM Rica",
            0.45: "Arquitetura Equilibrada",
            0.0: "Sono Eficiente e Compacto",
        },
        "chronotype": {
            0.65: "Noturno Genético",
            0.45: "Cronótipo Intermediário",
            0.0: "Matutino Genético",
        },
        "nocturnal": {
            0.65: "Explorador da Fronteira",
            0.45: "Sensibilidade Moderada",
            0.0: "Fronteira Definida",
        },
    }

    cat_profiles = profiles.get(category, profiles["dreamVividness"])
    for threshold in sorted(cat_profiles.keys(), reverse=True):
        if avg >= threshold:
            return cat_profiles[threshold]
    return "Típico"


# ---------------------------------------------------------------------------
# FUNÇÃO PRINCIPAL DE ANÁLISE
# ---------------------------------------------------------------------------

def analyze_dream_sleep(
    variants: Dict[str, Tuple[str, str, str]]
) -> Dict[str, Any]:
    """
    Executa a análise completa de Sonhos & Sono para as variantes do usuário.

    Args:
        variants: Dict mapeando rsid -> (chromosome, position, genotype)

    Returns:
        Dict com resultados completos da análise
    """
    dream_vividness = _analyze_snp_list(variants, DREAM_VIVIDNESS_SNPS)
    sleep_architecture = _analyze_snp_list(variants, SLEEP_ARCHITECTURE_SNPS)
    chronotype = _analyze_snp_list(variants, CHRONOTYPE_SNPS)
    nocturnal = _analyze_snp_list(variants, NOCTURNAL_PHENOMENA_SNPS)

    return {
        "dream_vividness": dream_vividness,
        "sleep_architecture": sleep_architecture,
        "chronotype": chronotype,
        "nocturnal_phenomena": nocturnal,
    }


def generate_dream_sleep_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Gera o JSON final para o relatório de Sonhos & Sono.

    Args:
        result: Output de analyze_dream_sleep()

    Returns:
        Dict compatível com o schema JSON do frontend GeneHealth
    """
    dv = result["dream_vividness"]
    sa = result["sleep_architecture"]
    ch = result["chronotype"]
    np_ = result["nocturnal_phenomena"]

    return {
        "reportType": "dream_sleep",
        "version": "1.0",
        "dreamVividness": {
            "summary": _generate_category_summary(dv, "dreamVividness"),
            "overallProfile": _generate_overall_profile(dv, "dreamVividness"),
            "traits": dv,
        },
        "sleepArchitecture": {
            "summary": _generate_category_summary(sa, "sleepArchitecture"),
            "overallProfile": _generate_overall_profile(sa, "sleepArchitecture"),
            "traits": sa,
        },
        "chronotype": {
            "summary": _generate_category_summary(ch, "chronotype"),
            "overallProfile": _generate_overall_profile(ch, "chronotype"),
            "traits": ch,
        },
        "nocturnalPhenomena": {
            "summary": _generate_category_summary(np_, "nocturnal"),
            "overallProfile": _generate_overall_profile(np_, "nocturnal"),
            "spiritualContext": SPIRITUAL_DREAM_CONTEXT,
            "traits": np_,
        },
        "disclaimers": DISCLAIMERS,
        "crossReferences": {
            "mindSpirit": (
                "HTR2A (rs6311, rs6313) e BDNF (rs6265) também aparecem no relatório "
                "Mente & Espírito. A mesma neurobiologia que amplifica a sensibilidade "
                "espiritual na vigília se manifesta como vida onírica intensa no sono — "
                "a consciência é um continuum."
            ),
        },
    }
