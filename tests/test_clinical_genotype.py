"""
Tests for the Live Clinical Data genotype annotation (scripts/clinical_genotype.py).

Run from the repo root (stdlib only, no pytest needed):

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m unittest discover -s tests -v

The central safety property: someone carrying two copies of the classified
(risk) allele must never be told they match the reference. The May 2026
version inferred "homozygous ⇒ most likely reference", which would have told
an HFE C282Y homozygote (rs1800562 AA) "Not a carrier signal".
"""

from __future__ import annotations

import io
import itertools
import json
import os
import sys
import unittest
from contextlib import redirect_stdout
from typing import Any, Dict, List, Optional
from unittest import mock

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from scripts import clinvar_live  # noqa: E402
from scripts.clinical_genotype import (  # noqa: E402
    annotate_user_genotypes,
    classify_genotype,
    parse_dbsnp_spdi,
)
from scripts.clinvar_live import enrich_with_clinvar  # noqa: E402
from scripts.run_full_analysis import parse_genome_file  # noqa: E402

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def z(genotype: str, ref: Optional[str], alts: List[str], variant: Optional[str] = None,
      block: Optional[str] = None) -> Dict[str, Any]:
    return classify_genotype(genotype, ref, alts, variant, reference_block=block)


# ── Pure classification ────────────────────────────────────────────────────


class HfeC282YTests(unittest.TestCase):
    """rs1800562: dbSNP ref G, alts A (C282Y, the ClinVar entry) and T."""

    REF, ALTS, VAR = "G", ["A", "T"], "A"

    def cls(self, g: str, block: Optional[str] = None) -> Dict[str, Any]:
        return z(g, self.REF, self.ALTS, self.VAR, block)

    def test_homozygous_risk_allele_is_never_reported_as_reference(self):
        r = self.cls("AA")
        self.assertEqual(r["zygosity"], "homozygous_alt")
        self.assertIs(r["likelyCarrier"], True)

    def test_homozygous_reference(self):
        r = self.cls("GG")
        self.assertEqual(r["zygosity"], "homozygous_ref")
        self.assertIs(r["likelyCarrier"], False)
        self.assertNotIn("strandFlipped", r)

    def test_heterozygous_either_order_and_separators(self):
        for g in ("GA", "AG", "G/A", "a g"):
            with self.subTest(g=g):
                r = self.cls(g)
                self.assertEqual(r["zygosity"], "heterozygous")
                self.assertIs(r["likelyCarrier"], True)

    def test_reverse_strand_file(self):
        self.assertEqual(self.cls("CC")["zygosity"], "homozygous_ref")
        self.assertTrue(self.cls("CC")["strandFlipped"])
        self.assertEqual(self.cls("CT")["zygosity"], "heterozygous")
        self.assertTrue(self.cls("CT")["strandFlipped"])

    def test_third_dbsnp_allele_is_not_reinterpreted(self):
        # TT is either a plus-strand C282F homozygote or a minus-strand C282Y
        # homozygote. Neither is reference; we refuse to guess.
        r = self.cls("TT")
        self.assertEqual(r["zygosity"], "unknown")
        self.assertEqual(r["zygosityReason"], "other_allele")
        self.assertIsNone(r["likelyCarrier"])

    def test_hemizygous_single_letter(self):
        self.assertEqual(self.cls("G")["zygosity"], "homozygous_ref")
        self.assertEqual(self.cls("A")["zygosity"], "homozygous_alt")

    def test_no_calls(self):
        for g in ("--", "-", "00", "NC", "NN", "", "  "):
            with self.subTest(g=g):
                r = self.cls(g)
                self.assertEqual(r["zygosity"], "no_call")
                self.assertIsNone(r["likelyCarrier"])

    def test_indel_and_garbage_genotypes(self):
        for g in ("II", "DI", "DD", "GGA", "AX"):
            with self.subTest(g=g):
                self.assertEqual(self.cls(g)["zygosity"], "unknown")

    def test_reference_block_only_affects_reference_calls(self):
        r = self.cls("GG", block="reference_allele_classified")
        self.assertEqual(r["zygosity"], "unknown")
        self.assertEqual(r["zygosityReason"], "reference_allele_classified")
        self.assertIsNone(r["likelyCarrier"])
        self.assertEqual(self.cls("GA", block="reference_allele_classified")["zygosity"], "heterozygous")
        self.assertEqual(self.cls("AA", block="reference_allele_classified")["zygosity"], "homozygous_alt")


class PalindromicTests(unittest.TestCase):
    def test_c_g_site_homozygotes_are_unknown(self):
        # rs1799945 (HFE H63D): ref C, ClinVar allele G. "CC" is reference on
        # the plus strand but H63D homozygous on the minus strand.
        for g in ("CC", "GG", "C", "G"):
            with self.subTest(g=g):
                r = z(g, "C", ["G", "T"], "G")
                self.assertEqual(r["zygosity"], "unknown")
                self.assertEqual(r["zygosityReason"], "strand_ambiguous")

    def test_c_g_site_heterozygote_is_still_heterozygous(self):
        for g in ("CG", "GC"):
            with self.subTest(g=g):
                self.assertEqual(z(g, "C", ["G", "T"], "G")["zygosity"], "heterozygous")

    def test_a_t_site(self):
        self.assertEqual(z("AA", "A", ["T"], "T")["zygosity"], "unknown")
        self.assertEqual(z("TT", "A", ["T"], "T")["zygosity"], "unknown")
        self.assertEqual(z("AT", "A", ["T"], "T")["zygosity"], "heterozygous")


class MultiAllelicTests(unittest.TestCase):
    def test_mthfr_like_four_allele_site_with_clinvar_allele(self):
        # rs1801133: dbSNP lists G>A,C,T; ClinVar classifies G>A (C677T).
        self.assertEqual(z("GG", "G", ["A", "C", "T"], "A")["zygosity"], "homozygous_ref")
        self.assertEqual(z("GA", "G", ["A", "C", "T"], "A")["zygosity"], "heterozygous")
        self.assertEqual(z("AA", "G", ["A", "C", "T"], "A")["zygosity"], "homozygous_alt")
        # CC could be a minus-strand GG or a plus-strand rare-C homozygote.
        self.assertEqual(z("CC", "G", ["A", "C", "T"], "A")["zygosity"], "unknown")
        self.assertEqual(z("TT", "G", ["A", "C", "T"], "A")["zygosity"], "unknown")

    def test_four_allele_site_without_a_known_variant_allele(self):
        # No ClinVar allele: GG could be plus GG (ref) or minus CC (alt).
        self.assertEqual(z("GG", "G", ["A", "C", "T"])["zygosity"], "unknown")
        # AA (plus) and TT (minus) are both non-reference homozygotes.
        self.assertEqual(z("AA", "G", ["A", "C", "T"])["zygosity"], "homozygous_alt")
        self.assertEqual(z("GA", "G", ["A", "C", "T"])["zygosity"], "unknown")

    def test_two_different_variant_alleles(self):
        r = z("AT", "G", ["A", "T"])
        self.assertEqual(r["zygosity"], "unknown")
        self.assertEqual(r["zygosityReason"], "two_different_variant_alleles")

    def test_missing_allele_data(self):
        self.assertEqual(z("GG", None, [])["zygosity"], "unknown")
        self.assertEqual(z("GG", "G", [])["zygosity"], "unknown")


class ExhaustiveSafetyTests(unittest.TestCase):
    """Brute-force every site configuration and genotype; check the invariant
    that matters: a genotype is called homozygous_ref only if, on the strand
    that reading used, every allele is the reference AND the other strand's
    reading cannot be a valid non-reference genotype."""

    GENOTYPES = ["".join(p) for p in itertools.product("ACGT", repeat=2)] + list("ACGT")

    def site_configs(self):
        bases = "ACGT"
        for ref in bases:
            others = [b for b in bases if b != ref]
            for n in (1, 2, 3):
                for alts in itertools.combinations(others, n):
                    for variant in list(alts) + [None]:
                        yield ref, list(alts), variant

    def test_variant_allele_never_yields_reference(self):
        checked = 0
        for ref, alts, variant in self.site_configs():
            site = {ref, *alts}
            for g in self.GENOTYPES:
                r = z(g, ref, alts, variant)
                checked += 1
                if r["zygosity"] != "homozygous_ref":
                    continue
                as_reported_ref = set(g) == {ref}
                as_flipped_ref = set(g) == {COMP[ref]}
                self.assertTrue(as_reported_ref or as_flipped_ref, (g, ref, alts, variant))
                if as_reported_ref:
                    # The opposite-strand reading (all comp(ref)) must not be
                    # a site allele that counts here, or it would be ambiguous.
                    counted = {ref, variant} if variant else site
                    self.assertNotIn(COMP[ref], counted, (g, ref, alts, variant))
                else:
                    # Flipped: the letters must not be a real allele at the site.
                    self.assertNotIn(COMP[ref], site, (g, ref, alts, variant))
                # The classified allele itself is never present.
                if variant:
                    self.assertNotIn(variant, g, (g, ref, alts, variant))
        self.assertGreater(checked, 1000)

    def test_two_copies_of_the_variant_on_the_reported_strand(self):
        for ref, alts, variant in self.site_configs():
            if not variant:
                continue
            r = z(variant * 2, ref, alts, variant)
            self.assertIn(r["zygosity"], ("homozygous_alt", "unknown"), (ref, alts, variant))


class ParseSpdiTests(unittest.TestCase):
    def test_multi_allelic(self):
        self.assertEqual(
            parse_dbsnp_spdi("NC_000006.12:26092912:G:A,NC_000006.12:26092912:G:T"),
            {"seq": "NC_000006.12", "pos": 26092912, "ref": "G", "alts": ["A", "T"]},
        )

    def test_rejects_indels_and_disagreement(self):
        self.assertIsNone(parse_dbsnp_spdi("NC_000001.11:100:GA:G"))
        self.assertIsNone(parse_dbsnp_spdi("NC_000001.11:100:G:A,NC_000001.11:101:G:T"))
        self.assertIsNone(parse_dbsnp_spdi(""))
        self.assertIsNone(parse_dbsnp_spdi(None))


# ── Pipeline: fake NCBI → enrich_with_clinvar → annotate ───────────────────


def _dbsnp(rs_digits: str, gene: str, spdi: str, sig: str) -> Dict[str, Any]:
    return {
        "SNP_ID": rs_digits,
        "CLINICAL_SIGNIFICANCE": sig,
        "GENES": [{"NAME": gene}],
        "SNP_CLASS": "snv",
        "FXN_CLASS": "missense_variant",
        "SPDI": spdi,
    }


def _clinvar(uid: str, title: str, desc: str, entries: List[tuple], traits=("Some condition",),
             review: str = "criteria provided, multiple submitters, no conflicts") -> Dict[str, Any]:
    return {
        "uid": uid,
        "title": title,
        "germline_classification": {
            "description": desc,
            "review_status": review,
            "last_evaluated": "2025/01/02 00:00",
            "trait_set": [{"trait_name": t} for t in traits],
        },
        "variation_set": [
            {
                "canonical_spdi": spdi,
                "variation_xrefs": [{"db_source": "dbSNP", "db_id": rs}] if rs else [],
            }
            for rs, spdi in entries
        ],
    }


class FakeNCBI:
    """Stands in for NCBIClient with synthetic dbSNP + ClinVar records."""

    DBSNP = {
        # HFE C282Y — multi-allelic, simple ClinVar SNV record.
        "1800562": _dbsnp("1800562", "HFE", "NC_000006.12:26092912:G:A,NC_000006.12:26092912:G:T", "pathogenic,risk-factor"),
        # HFE H63D — C/G palindromic with the ClinVar allele.
        "1799945": _dbsnp("1799945", "HFE", "NC_000006.12:26090950:C:G,NC_000006.12:26090950:C:T", "pathogenic"),
        # COMT — simple biallelic drug-response site.
        "4680": _dbsnp("4680", "COMT", "NC_000022.11:19963747:G:A", "benign,drug-response"),
        # APOE — haplotype record + a non-benign reference-equal record.
        "429358": _dbsnp("429358", "APOE", "NC_000019.10:44908683:T:C", "pathogenic,risk-factor"),
        # Sites whose ClinVar data disagrees with dbSNP / is missing.
        "1000001": _dbsnp("1000001", "GENE1", "NC_000001.11:500:G:A", "pathogenic"),
        "1000002": _dbsnp("1000002", "GENE2", "NC_000001.11:600:C:T", "pathogenic"),
        # Informational site with no ClinVar call (dbSNP single ALT).
        "1000003": _dbsnp("1000003", "GENE3", "NC_000001.11:700:A:G", "benign"),
        # Indel site — genotype letters cannot be compared.
        "1000004": _dbsnp("1000004", "GENE4", "NC_000001.11:800:GA:G", "pathogenic"),
        # Review findings (2026-11-23):
        "2000001": _dbsnp("2000001", "REF1", "NC_000002.12:10:T:C", "pathogenic"),
        "2000002": _dbsnp("2000002", "REF2", "NC_000002.12:20:T:C", "pathogenic"),
        "2000003": _dbsnp("2000003", "CYP", "NC_000010.11:430:C:T", "drug-response"),
        "2000004": _dbsnp("2000004", "MANY", "NC_000002.12:40:G:A", "pathogenic"),
        "2000005": _dbsnp("2000005", "LOWPEN", "NC_000002.12:50:G:A", "pathogenic-low-penetrance"),
        "2000006": _dbsnp("2000006", "NOCLASS", "NC_000002.12:60:G:A", "pathogenic"),
        "2000007": _dbsnp("2000007", "HAPREF", "NC_000002.12:70:T:C", "pathogenic"),
        "2000008": _dbsnp("2000008", "OTHER", "NC_000002.12:80:G:A", "other"),
    }

    CLINVAR_SEARCH = {
        "rs1800562": ["9"],
        "rs1799945": ["10"],
        "rs4680": ["17591"],
        "rs429358": ["441268", "17863"],
        "rs1000001": ["900001"],
        "rs1000002": [],
        "rs1000004": [],
        "rs2000001": ["910001", "910002"],
        "rs2000002": ["920001"],
        "rs2000003": ["930001", "930002"],
        # A full page of hits: not every record could be read.
        "rs2000004": ["940001"] + [str(950000 + i) for i in range(49)],
        "rs2000006": ["960001", "960002"],
        "rs2000007": ["970001"],
    }

    CLINVAR = {
        "9": _clinvar("9", "NM_000410.4(HFE):c.845G>A (p.Cys282Tyr)", "Pathogenic",
                      [("1800562", "NC_000006.12:26092912:G:A")], ("Hereditary hemochromatosis",)),
        "10": _clinvar("10", "NM_000410.4(HFE):c.187C>G (p.His63Asp)", "Pathogenic",
                       [("1799945", "NC_000006.12:26090950:C:G")]),
        "17591": _clinvar("17591", "NM_000754.4(COMT):c.472G>A (p.Val158Met)", "drug response",
                          [("4680", "NC_000022.11:19963747:G:A")]),
        # Haplotype: two entries; ours is the rs429358 one.
        "441268": _clinvar("441268", "NM_000041.3(APOE):c.[137T>C;388T>C]", "Pathogenic/Likely pathogenic",
                           [("769452", "NC_000019.10:44907852:T:C"), ("429358", "NC_000019.10:44908683:T:C")]),
        # Reference-equal record with a non-benign classification.
        "17863": _clinvar("17863", "NM_000041.4(APOE):c.[388=;526=]", "Pathogenic",
                          [("429358", "NC_000019.10:44908683:T:T")], ("APOE3 ISOFORM",)),
        # ClinVar SPDI at a different position than dbSNP → inconsistent.
        "900001": _clinvar("900001", "NM_X:c.1G>A", "Pathogenic", [("1000001", "NC_000001.11:501:G:A")]),
        # A reference-equal record that outranks the real change (more stars,
        # more conditions) must not be the one shown on the card.
        "910001": _clinvar("910001", "NM_R1:c.10=", "Pathogenic", [("2000001", "NC_000002.12:10:T:T")],
                           ("Ref condition A", "Ref condition B", "Ref condition C"),
                           review="reviewed by expert panel"),
        "910002": _clinvar("910002", "NM_R1:c.10T>C", "Pathogenic", [("2000001", "NC_000002.12:10:T:C")],
                           ("Alt condition",), review="criteria provided, single submitter"),
        # Only a reference-equal record exists.
        "920001": _clinvar("920001", "NM_R2:c.20=", "Pathogenic", [("2000002", "NC_000002.12:20:T:T")]),
        # Haplotype whose real change is elsewhere: the reference base at our
        # site is context, not what is classified.
        "930001": _clinvar("930001", "NM_CYP:c.[430C=;1075A>C]", "drug response",
                           [("2000003", "NC_000010.11:430:C:C"), ("2999999", "NC_000010.11:1075:A:C")]),
        "930002": _clinvar("930002", "NM_CYP:c.430C>T", "drug response", [("2000003", "NC_000010.11:430:C:T")]),
        "940001": _clinvar("940001", "NM_MANY:c.40G>A", "Pathogenic", [("2000004", "NC_000002.12:40:G:A")]),
        # Reference-equal record with no classification: not a reason to doubt.
        "960001": _clinvar("960001", "NM_NC:c.60=", "no classification for the single variant",
                           [("2000006", "NC_000002.12:60:G:G")]),
        "960002": _clinvar("960002", "NM_NC:c.60G>A", "Pathogenic", [("2000006", "NC_000002.12:60:G:A")]),
        # The only record is a pathogenic haplotype carrying the reference
        # base at our site: it is the one shown, so its classification is
        # about the allele the user would have two copies of.
        "970001": _clinvar("970001", "NM_HR:c.[70=;500G>A]", "Pathogenic",
                           [("2000007", "NC_000002.12:70:T:T"), ("2777777", "NC_000002.12:500:G:A")]),
    }

    def esummary(self, db: str, ids):
        src = self.DBSNP if db == "snp" else self.CLINVAR
        return [dict(src[i]) for i in ids if i in src]

    def esearch(self, db: str, term: str, retmax: int = 10):
        return list(self.CLINVAR_SEARCH.get(term.rstrip("*"), []))[:retmax]


GENOME_23ANDME = "\n".join([
    "# This data file generated by 23andMe",
    "# rsid\tchromosome\tposition\tgenotype",
    "rs1800562\t6\t26093141\tAA",   # C282Y homozygote (the original bug)
    "rs1799945\t6\t26091179\tCC",   # palindromic homozygote
    "rs4680\t22\t19951271\tAG",     # heterozygote
    "rs429358\t19\t45411941\tTT",   # reference, but ClinVar classifies the reference allele
    "rs1000001\t1\t500\tGG",        # NCBI sources disagree
    "rs1000002\t1\t600\tCC",        # high severity, no ClinVar record found
    "rs1000003\t1\t700\tAA",        # informational, homozygous reference
    "rs1000004\t1\t800\tGG",        # indel site
    "rs9999999\t1\t900\tAG",        # not enriched
])

GENOME_ANCESTRY = "\r\n".join([
    "#AncestryDNA raw data download",
    "rsid\tchromosome\tposition\tallele1\tallele2",
    "rs1800562\t6\t26093141\tG\tG",
    "rs1799945\t6\t26091179\tC\tG",
    "rs4680\t22\t19951271\t0\t0",
    "rs429358\t19\t45411941\tC\tC",
    "rs1000003\t1\t700\tA\tG",
])

ALL_RSIDS = ["rs1800562", "rs1799945", "rs4680", "rs429358", "rs1000001", "rs1000002", "rs1000003", "rs1000004"]


class PipelineTests(unittest.TestCase):
    def enrich(self) -> Dict[str, Dict[str, Any]]:
        return enrich_with_clinvar(ALL_RSIDS, client=FakeNCBI(), max_rsids=40)

    def test_enrichment_alleles(self):
        cl = self.enrich()
        self.assertEqual(cl["rs1800562"]["refAllele"], "G")
        self.assertEqual(cl["rs1800562"]["altAlleles"], ["A", "T"])
        self.assertEqual(cl["rs1800562"]["variantAllele"], "A")
        # Haplotype record: the allele comes from the entry for *this* rsID.
        self.assertEqual(cl["rs429358"]["variantAllele"], "C")
        self.assertTrue(cl["rs429358"]["referenceAlleleClassified"])
        # Disagreeing NCBI sources: allele fields dropped.
        for key in ("refAllele", "altAlleles", "variantAllele"):
            self.assertNotIn(key, cl["rs1000001"])
        # Informational, single dbSNP ALT used as the variant allele.
        self.assertEqual(cl["rs1000003"]["variantAllele"], "G")
        # Indel site: no allele fields at all.
        self.assertNotIn("refAllele", cl["rs1000004"])
        # Internal bookkeeping never leaks into the report JSON.
        for payload in cl.values():
            self.assertFalse([k for k in payload if k.startswith("_")], payload)

    def test_23andme_genome(self):
        cl = self.enrich()
        annotate_user_genotypes(cl, parse_genome_file(GENOME_23ANDME, "23andme"))
        got = {rs: (p["userGenotype"], p["zygosity"], p.get("zygosityReason")) for rs, p in cl.items()}
        self.assertEqual(got, {
            "rs1800562": ("AA", "homozygous_alt", None),
            "rs1799945": ("CC", "unknown", "strand_ambiguous"),
            "rs4680": ("AG", "heterozygous", None),
            "rs429358": ("TT", "unknown", "reference_allele_classified"),
            "rs1000001": ("GG", "unknown", "no_allele_data"),
            "rs1000002": ("CC", "unknown", "reference_not_confirmed"),
            "rs1000003": ("AA", "homozygous_ref", None),
            "rs1000004": ("GG", "unknown", "no_allele_data"),
        })
        self.assertIs(cl["rs1800562"]["likelyCarrier"], True)
        self.assertIs(cl["rs1000003"]["likelyCarrier"], False)
        self.assertIsNone(cl["rs1799945"]["likelyCarrier"])
        json.dumps(cl)  # must stay JSON-serialisable

    def test_ancestry_genome_with_crlf_and_missing_rows(self):
        cl = self.enrich()
        annotate_user_genotypes(cl, parse_genome_file(GENOME_ANCESTRY, "ancestry"))
        self.assertEqual(cl["rs1800562"]["zygosity"], "homozygous_ref")
        self.assertEqual(cl["rs1799945"]["zygosity"], "heterozygous")
        self.assertEqual(cl["rs4680"]["zygosity"], "no_call")
        self.assertEqual(cl["rs4680"]["userGenotype"], "--")
        self.assertEqual(cl["rs429358"]["zygosity"], "homozygous_alt")
        self.assertEqual(cl["rs1000003"]["zygosity"], "heterozygous")
        # Not in the file at all.
        self.assertEqual(cl["rs1000002"]["zygosity"], "no_call")
        self.assertEqual(cl["rs1000002"]["zygosityReason"], "not_in_file")
        self.assertIsNone(cl["rs1000002"]["userGenotype"])

    def test_one_bad_card_does_not_break_the_rest(self):
        cl = self.enrich()
        cl["rs4680"]["altAlleles"] = 12345  # corrupt one card
        annotate_user_genotypes(cl, parse_genome_file(GENOME_23ANDME, "23andme"))
        self.assertEqual(cl["rs4680"]["zygosity"], "unknown")
        self.assertEqual(cl["rs1800562"]["zygosity"], "homozygous_alt")


class ReviewFindingTests(unittest.TestCase):
    """Regression tests for the 2026-11-23 adversarial review."""

    RSIDS = ["rs2000001", "rs2000002", "rs2000003", "rs2000004", "rs2000005", "rs2000006",
             "rs2000007", "rs2000008"]

    def annotated(self, genotypes: Dict[str, str]) -> Dict[str, Dict[str, Any]]:
        cl = enrich_with_clinvar(self.RSIDS, client=FakeNCBI(), max_rsids=40)
        annotate_user_genotypes(cl, {rs: ("1", "1", g) for rs, g in genotypes.items()})
        return cl

    def test_reference_record_is_not_shown_when_a_real_change_exists(self):
        cl = self.annotated({"rs2000001": "TT"})
        card = cl["rs2000001"]
        self.assertEqual(card["clinvarTitle"], "NM_R1:c.10T>C")
        self.assertEqual(card["variantAllele"], "C")
        self.assertTrue(card["referenceAlleleClassified"])
        self.assertEqual(card["zygosity"], "unknown")
        self.assertEqual(card["zygosityReason"], "reference_allele_classified")

    def test_only_a_reference_record_drops_the_variant_allele(self):
        cl = self.annotated({"rs2000002": "TT"})
        card = cl["rs2000002"]
        self.assertEqual(card["clinvarTitle"], "NM_R2:c.20=")
        self.assertNotIn("variantAllele", card)
        self.assertEqual(card["zygosity"], "unknown")
        self.assertEqual(card["zygosityReason"], "reference_allele_classified")
        cl = self.annotated({"rs2000002": "TC"})
        self.assertEqual(cl["rs2000002"]["zygosity"], "heterozygous")

    def test_haplotype_context_does_not_block_a_reference_call(self):
        cl = self.annotated({"rs2000003": "CC"})
        card = cl["rs2000003"]
        self.assertNotIn("referenceAlleleClassified", card)
        self.assertEqual(card["variantAllele"], "T")
        self.assertEqual(card["zygosity"], "homozygous_ref")

    def test_truncated_clinvar_search_is_not_a_confirmation(self):
        cl = self.annotated({"rs2000004": "GG"})
        card = cl["rs2000004"]
        self.assertTrue(card["referenceCheckIncomplete"])
        self.assertEqual(card["zygosity"], "unknown")
        self.assertEqual(card["zygosityReason"], "reference_not_confirmed")
        self.assertEqual(self.annotated({"rs2000004": "GA"})["rs2000004"]["zygosity"], "heterozygous")

    def test_risk_label_outside_the_tier_table_needs_clinvar(self):
        cl = self.annotated({"rs2000005": "GG"})
        card = cl["rs2000005"]
        self.assertEqual(card["severity"], "info")  # tier table does not know the label
        self.assertNotIn("clinvarId", card)
        self.assertEqual(card["zygosity"], "unknown")
        self.assertEqual(card["zygosityReason"], "reference_not_confirmed")

    def test_displayed_haplotype_with_reference_base_blocks(self):
        cl = self.annotated({"rs2000007": "TT"})
        card = cl["rs2000007"]
        self.assertEqual(card["clinvarTitle"], "NM_HR:c.[70=;500G>A]")
        self.assertNotIn("variantAllele", card)
        self.assertTrue(card["referenceAlleleClassified"])
        self.assertEqual(card["zygosity"], "unknown")
        self.assertEqual(card["zygosityReason"], "reference_allele_classified")

    def test_other_label_without_clinvar_is_not_an_all_clear(self):
        card = self.annotated({"rs2000008": "GG"})["rs2000008"]
        self.assertEqual(card["zygosity"], "unknown")
        self.assertEqual(card["zygosityReason"], "reference_not_confirmed")

    def test_unclassified_reference_record_does_not_block(self):
        cl = self.annotated({"rs2000006": "GG"})
        card = cl["rs2000006"]
        self.assertNotIn("referenceAlleleClassified", card)
        self.assertEqual(card["zygosity"], "homozygous_ref")


class RunAnalysisEndToEndTests(unittest.TestCase):
    """Full run_analysis() on a synthetic 23andMe genome with NCBI faked."""

    def test_clinical_live_report_carries_genotypes(self):
        from scripts import run_full_analysis

        with mock.patch.object(clinvar_live, "default_client", return_value=FakeNCBI()), \
                mock.patch.dict(os.environ, {"ENABLE_CLINVAR_LIVE": "true"}), \
                redirect_stdout(io.StringIO()):
            result = run_full_analysis.run_analysis(GENOME_23ANDME, "23andme")

        report = json.loads(result["reports"]["clinical_live"])
        variants = report["variants"]
        self.assertIn("rs1800562", variants, sorted(variants))
        hfe = variants["rs1800562"]
        self.assertEqual(hfe["userGenotype"], "AA")
        self.assertEqual(hfe["zygosity"], "homozygous_alt")
        self.assertIs(hfe["likelyCarrier"], True)
        self.assertEqual(hfe["variantAllele"], "A")
        for rs, p in variants.items():
            self.assertIn(p["zygosity"], ("homozygous_ref", "heterozygous", "homozygous_alt", "no_call", "unknown"), rs)


if __name__ == "__main__":
    unittest.main()
