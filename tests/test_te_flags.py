import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from te_flags import te_class  # noqa: E402

TE = [
    "Plavaka: Plavaka transposase", "gag-asp_proteas: gag-polyprotein putative aspartyl protease",
    "RVP_2: Retroviral aspartyl protease", "Peptidase_A2B: Ty3 transposon peptidase",
    "Ty3_capsid: Ty3 transposon capsid-like protein", "PRO_PEG10-RTL1: PEG10/RTL1 protease domain",
    "PEG10_N-capsid: Retrotransposon-derived protein PEG10, N-terminal capsid-like domain",
    "MULE: MULE transposase domain", "Dimer_Tnp_hAT: hAT family C-terminal dimerisation region",
    "HTH_OrfB_IS605: Helix-turn-helix domain", "OrfB_IS605: Probable transposase",
    "TNP-like_RNaseH_N: TNP-like, RNase H N-terminal domain", "Cas12f1-like_TNB: Cas12f1-like, TNB domain",
    "RVT_1: Reverse transcriptase (RNA-dependent DNA polymerase)", "RVT_2: Reverse transcriptase (RNA-dependent DNA polymerase)",
    "RT_RNaseH_2: RNase H-like domain found in reverse transcriptase", "zf-RVT: zinc-binding in reverse transcriptase",
    "Exo_endo_phos_2: Endonuclease-reverse transcriptase", "rve: Integrase core domain", "rva_4: Integrase core domain",
    "Integrase_H2C2: Integrase zinc binding domain", "gag_pre-integrs: GAG-pre-integrase domain",
    "Retrotran_gag_2: gag-polypeptide of LTR copia-type", "Helitron_like_N: Helitron helicase-like domain at N-terminus",
    "KDZ: Kyakuja-Dileera-Zisupton transposase", "CxC1: CxC1 like cysteine cluster associated with KDZ transposases",
    "DDE_1: DDE superfamily endonuclease", "DDE_Tnp_4: DDE superfamily endonuclease",
    "HTH_Tnp_Tc5: Tc5 transposase DNA-binding domain", "Pol_BBD: Pol polyprotein, beta-barrel domain",
    "SH3_retrovirus: Retroviral polymerase SH3-like domain", "Spuma_A9PTase: Spumavirus aspartic protease (A9)",
    "Plant_tran: Plant transposon protein", "HTH_Tnp_Tc3_2: Transposase", "Transposase_21: Transposase family tnp2",
]
OTHER = [
    "Zn_clus: Fungal Zn(2)-Cys(6) binuclear cluster domain", "F-box: F-box domain", "F-box-like: F-box-like",
    "Vhr1: Transcription factor Vhr1", "RNase_H: RNase H", "Asp_protease_2: Aspartyl protease",
    "Asp_protease: Aspartyl protease", "Uma2: Putative restriction endonuclease",
    "Mrr_cat: Restriction endonuclease", "LAGLIDADG_1: LAGLIDADG endonuclease", "FAR1: FAR1 DNA-binding domain",
    "Globin: Globin", "Fasciclin: Fasciclin domain", "BTB: BTB/POZ domain", "Glyco_hydro_16: Glycosyl hydrolases family 16",
    "DNA_pol_B: DNA polymerase family B", "CFEM: CFEM domain", "HMG_box: HMG (high mobility group) box",
]
REPEAT = ["Ank_2: Ankyrin repeats (3 copies)", "WD40: WD domain, G-beta repeat", "Sel1: Sel1 repeat",
          "LRR_8: Leucine rich repeat", "zf-CCHC: Zinc knuckle", "Exo_endo_phos: Endonuclease/Exonuclease/phosphatase family"]


@pytest.mark.parametrize("name", TE)
def test_te(name):
    assert te_class(name) == "TE"


@pytest.mark.parametrize("name", OTHER)
def test_other(name):
    assert te_class(name) == "other"


@pytest.mark.parametrize("name", REPEAT)
def test_repeat(name):
    assert te_class(name) == "repeat"


def test_name_without_description():
    assert te_class("Plavaka") == "TE"
    assert te_class("") == "other"
