"""Coordinator: identify drug -> retrieve from each source -> normalize -> compare -> flag gaps/conflicts."""
import pandas as pd
from models import DrugIdentity, DrugReport
import sources as S

PK_PD_NOTE = ("No PK studies for this drug in PK-DB, and PubChem/PharmGKB/Open Targets carry no PK parameters. "
              "PK/PD numbers in the simulator come from the curated library (pk_data.py).")
PKDB_NOTE = ("PK-DB lists studies for this drug, but parameter values (AUC, clearance, half-life) are not returned "
             "to anonymous API users; values used by the simulator come from the curated library.")


def run(drug: str, source_fns=None) -> DrugReport:
    identity = DrugIdentity(query=drug.strip())
    report = DrugReport(identity=identity)
    for fn in (source_fns or S.SOURCES):
        ev, status = fn(identity)
        report.evidence += ev
        report.source_status.append(status)

    pgx = [e for e in report.evidence if e.kind in ("pgx_annotation", "gene_drug")]
    report.genes_of_interest = sorted({g for e in report.evidence if e.gene for g in e.gene.split(",")})
    _analyze(report, pgx)
    return report


def _analyze(report: DrugReport, pgx):
    # Gaps
    for s in report.source_status:
        if not s.ok:
            report.gaps.append(f"{s.source} failed: {s.error}")
        elif s.n_records == 0:
            report.gaps.append(f"{s.source} returned no records ({s.error or 'no match'})")
    report.gaps.append(PKDB_NOTE if any(e.kind == "pk_study" for e in report.evidence) else PK_PD_NOTE)
    if not any(e.population for e in pgx):
        report.gaps.append("No ancestry/population information retrieved for any gene-drug finding. "
                           "Note: genetic ancestry and ethnicity are not interchangeable.")
    if not pgx:
        report.gaps.append("No pharmacogenomic variant evidence found for this drug.")

    # Conflicts / agreement: genes flagged by one source but not another
    by_src = {}
    for e in pgx:
        for g in (e.gene or "").split(","):
            if g:
                by_src.setdefault(e.source, set()).add(g)
    if len(by_src) >= 2:
        allg = set().union(*by_src.values())
        for g in sorted(allg):
            have = [s for s, gs in by_src.items() if g in gs]
            if len(have) < len(by_src):
                report.conflicts.append(f"{g}: reported by {', '.join(have)} only")


def evidence_df(report: DrugReport) -> pd.DataFrame:
    return pd.DataFrame([e.model_dump() for e in report.evidence])


def to_markdown(report: DrugReport) -> str:
    i = report.identity
    lines = [f"# Drug response report: {i.name or i.query}", "",
             f"- PubChem CID: {i.pubchem_cid} | ChEMBL: {i.chembl_id} | PharmGKB: {i.pharmgkb_id}",
             f"- Formula: {i.formula} | MW: {i.molecular_weight}", "",
             "## Source status"]
    lines += [f"- {s.source}: {'ok' if s.ok else 'FAILED'} ({s.n_records} records){' - ' + s.error if s.error else ''}"
              for s in report.source_status]
    pk = [e for e in report.evidence if e.kind == "pk_study"]
    lines += ["", "## Q1. Absorption / distribution / metabolism / elimination"]
    lines += [f"- [PK-DB] {e.summary} {e.url}" for e in pk] or [PK_PD_NOTE]
    lines += [
              "", "## Q2. Genetic variants affecting metabolism or response"]
    pgx = [e for e in report.evidence if e.kind == "pgx_annotation"]
    lines += [f"- [{e.source}] {e.gene or '?'} {e.variant or ''} (level {e.evidence_level or 'n/a'}): {e.summary}"
              for e in pgx] or ["- none retrieved"]
    lines += ["", "## Mechanism / targets"]
    lines += [f"- [{e.source}] {e.gene or '?'}: {e.summary}" for e in report.evidence if e.kind == "mechanism"] or ["- none"]
    lines += ["", "## Q3. Demographic / ancestry evidence"]
    pops = [e for e in pgx if e.population]
    lines += [f"- {e.gene}: {e.population}" for e in pops] or ["- not available from current sources"]
    lines += ["", "## Conflicts across sources"] + ([f"- {c}" for c in report.conflicts] or ["- none detected"])
    lines += ["", "## What remains uncertain / missing"] + [f"- {g}" for g in report.gaps]
    return "\n".join(lines)


def pkpd_markdown(profile, sim, flag_list, report: DrugReport | None = None) -> str:
    """Full drug-response report: patient-specific PK/PD + (optional) live evidence."""
    from pk_data import DRUGS
    d = DRUGS[profile.drug]
    m, mr, u = sim.metrics, sim.metrics_ref, sim.unit
    per_day = profile.doses_per_day or round(24 / profile.interval_h)
    L = [f"# Drug response report: {profile.drug} ({d['brand']})", "",
         "_Research prototype. Not for clinical use._", "",
         "## Patient profile",
         f"- {profile.age} y, {profile.sex}, {profile.weight_kg:.0f} kg; eGFR {profile.egfr:.0f}; liver: {profile.hepatic}",
         f"- Regimen: {profile.dose_mg:g} mg x {per_day}/day (every {profile.interval_h:g} h), adherence {profile.adherence:.0%}",
         f"- Genotype: {', '.join(f'{g} {p}' for g, p in profile.phenotypes.items()) or 'not provided'}",
         f"- Co-medications / supplements: {', '.join(profile.comeds) or 'none'}", "",
         "## Q1. Absorption, distribution, metabolism, elimination",
         f"- Typical: t1/2 {d['t_half']:g} h, Tmax {d['tmax']:g} h, V/F {d['vd_f']:g} L/kg, protein binding {d['protein_binding']}, "
         f"renal fraction {d['fe']:.0%}; pathways: {', '.join(f'{k} {v:.0%}' for k, v in d['fm'].items()) or 'renal'}",
         f"- This patient: t1/2 {sim.t_half_patient:.1f} h; clearance {sim.clearance.cl_patient:.2f} L/h "
         f"(typical {sim.clearance.cl_ref:.2f}); exposure ~{sim.clearance.ratio:.2f}x typical",
         f"- Steady state (last simulated day): Cmax {m['cmax_ss']:.3g}, Cmin {m['cmin_ss']:.3g}, Cavg {m['cavg_ss']:.3g} {u} "
         f"(typical adult Cavg {mr['cavg_ss']:.3g} {u})"]
    if d.get("range"):
        L.append(f"- Reference range {d['range'][0]}-{d['range'][1]} {u}; time in range over last week: {m.get('pct_in_range', 0):.0f}%")
    L += ["", "## Q2. Genetic variants affecting metabolism or response"]
    L += [f"- {g}: {profile.phenotypes.get(g, 'not provided')}" for g in d["pgx_genes"]] or ["- No actionable PK genes for this drug"]
    L += [f"- PD genetics: {d['pd_genes']}"]
    if report:
        pgx = [e for e in report.evidence if e.kind == "pgx_annotation" and e.source == "PharmGKB"]
        L += [f"- [ClinPGx {e.evidence_level or 'n/a'}] {e.summary}" for e in pgx[:15]]
    L += ["", "## Q3. Demographic / ancestry evidence",
          "- Allele frequencies differ by ancestry (e.g. CYP2C19 no-function alleles are more common in East Asian populations; "
          "CYP2D6 UM more common in some North African / Middle Eastern populations).",
          "- Genetic ancestry and ethnicity are not interchangeable."]
    if report:
        pops = [e for e in report.evidence if e.population]
        L += [f"- {e.gene}: {e.population}" for e in pops] or ["- No population-level detail retrieved from live sources."]
    L += ["", "## Q4. Factors that explain different responses / duration of action"]
    L += [f"- {lbl}: clearance x{f:.2f}" for lbl, f in sim.clearance.factors] or ["- None beyond typical adult"]
    L += ["", "## Pharmacodynamics", f"- Mechanism: {d['mechanism']}", f"- Onset: {d['onset']}",
          "- Adverse effects to watch: " + "; ".join(d["adverse"])]
    L += ["", "## Flags"] + [f"- **{f.level.upper()}** [{f.area}] {f.title}: {f.detail}" for f in flag_list]
    if report:
        L += ["", "## Conflicts across sources"] + ([f"- {c}" for c in report.conflicts[:20]] or ["- none detected"])
        L += ["", "## What remains uncertain / missing"] + [f"- {g}" for g in report.gaps]
    L += ["", f"_Sources for PK values: {d['source']}_"]
    return "\n".join(L)
