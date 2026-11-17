"""Shared data models: every source gets normalized into these."""
from typing import Optional
from pydantic import BaseModel, Field


class DrugIdentity(BaseModel):
    query: str
    name: Optional[str] = None
    pubchem_cid: Optional[int] = None
    chembl_id: Optional[str] = None
    pharmgkb_id: Optional[str] = None
    formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    smiles: Optional[str] = None
    inchikey: Optional[str] = None


class Evidence(BaseModel):
    """One normalized finding from one source (a gene-drug association, a mechanism, etc.)."""
    source: str                      # "PubChem" | "PharmGKB" | "OpenTargets"
    kind: str                        # "property" | "gene_drug" | "mechanism" | "pgx_annotation"
    gene: Optional[str] = None
    variant: Optional[str] = None    # rsID or star allele
    summary: str
    evidence_level: Optional[str] = None   # e.g. PharmGKB level 1A, 2B...
    population: Optional[str] = None       # ancestry/ethnicity if the source reports it
    url: Optional[str] = None


class SourceStatus(BaseModel):
    source: str
    ok: bool
    n_records: int = 0
    error: Optional[str] = None


class DrugReport(BaseModel):
    identity: DrugIdentity
    evidence: list[Evidence] = Field(default_factory=list)
    source_status: list[SourceStatus] = Field(default_factory=list)
    genes_of_interest: list[str] = Field(default_factory=list)
    conflicts: list[str] = Field(default_factory=list)
    gaps: list[str] = Field(default_factory=list)
