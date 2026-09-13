"""Shared data structures for the pangenome cluster-profile analysis
(notes/superpowers/specs/2026-09-13-pangenome-cluster-profile-design.md).

PresenceMatrix holds one row per gene family (a tier-1 mmseqs/diamond
cluster representative) and one column per strain, with a three-state
call per cell -- see the Global Constraints in this plan's header for
what each state means.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

PRESENT = "present"
GENOME_ONLY = "genome_only"
ABSENT = "absent"
STATES = (PRESENT, GENOME_ONLY, ABSENT)


def read_cluster_tsv(path: str | Path) -> dict[str, str]:
    """Parse an mmseqs/diamond cluster TSV (rep\tmember per line) into
    {member_id: rep_id}."""
    member_to_rep: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            member_to_rep[member] = rep
    return member_to_rep


def build_families(member_to_rep: dict[str, str]) -> dict[str, list[str]]:
    """Return {rep_id: sorted [member_ids]} for every cluster, including
    singletons -- unlike NovInvenio's lib/clusters.py::build_families,
    which drops singleton clusters (irrelevant there; here a singleton
    family is a real, reportable frequency bin)."""
    families: dict[str, list[str]] = {}
    for member, rep in member_to_rep.items():
        families.setdefault(rep, []).append(member)
    return {rep: sorted(members) for rep, members in families.items()}


@dataclass
class PresenceMatrix:
    families: list[str]
    strains: list[str]
    calls: dict[tuple[str, str], str] = field(default_factory=dict)
    copy_number: dict[tuple[str, str], int] = field(default_factory=dict)

    def set_call(self, family: str, strain: str, state: str, copies: int = 0) -> None:
        if state not in STATES:
            raise ValueError(f"invalid state {state!r}, must be one of {STATES}")
        self.calls[(family, strain)] = state
        if copies:
            self.copy_number[(family, strain)] = copies

    def call(self, family: str, strain: str) -> str:
        return self.calls.get((family, strain), ABSENT)

    def is_present(self, family: str, strain: str) -> bool:
        return self.call(family, strain) in (PRESENT, GENOME_ONLY)

    def presence_vector(self, family: str) -> list[bool]:
        return [self.is_present(family, s) for s in self.strains]

    def strain_count(self, family: str) -> int:
        return sum(self.presence_vector(family))

    def frequency(self, family: str) -> float:
        if not self.strains:
            return 0.0
        return self.strain_count(family) / len(self.strains)

    def to_tsv(self, path: str | Path) -> None:
        with open(path, "w") as fh:
            fh.write("family\t" + "\t".join(self.strains) + "\n")
            for fam in self.families:
                row = [self.call(fam, s) for s in self.strains]
                fh.write(fam + "\t" + "\t".join(row) + "\n")

    @classmethod
    def from_tsv(cls, path: str | Path) -> "PresenceMatrix":
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            strains = header[1:]
            families: list[str] = []
            calls: dict[tuple[str, str], str] = {}
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                fam = parts[0]
                families.append(fam)
                for strain, state in zip(strains, parts[1:]):
                    calls[(fam, strain)] = state
        pm = cls(families=families, strains=strains)
        pm.calls = calls
        return pm
