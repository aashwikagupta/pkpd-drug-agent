"""API clients. Each returns normalized Evidence and never raises: failures are reported as SourceStatus.

NOTE: written without live API access. Endpoint shapes are from memory of the public docs;
verify each against the docs and expect to fix field names on first real run.
"""
import requests
from models import DrugIdentity, Evidence, SourceStatus

TIMEOUT = 20
HEADERS = {"User-Agent": "pgx-agent-prototype/0.1"}


# ---------- PubChem (identity + chemical properties) ----------
def pubchem(identity: DrugIdentity) -> tuple[list[Evidence], SourceStatus]:
    try:
        props = "MolecularFormula,MolecularWeight,SMILES,InChIKey"
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{identity.query}/property/{props}/JSON"
        r = requests.get(url, timeout=TIMEOUT, headers=HEADERS)
        r.raise_for_status()
        p = r.json()["PropertyTable"]["Properties"][0]
        identity.pubchem_cid = p.get("CID")
        identity.formula = p.get("MolecularFormula")
        identity.molecular_weight = float(p["MolecularWeight"]) if p.get("MolecularWeight") else None
        identity.smiles = p.get("SMILES") or p.get("CanonicalSMILES") or p.get("ConnectivitySMILES")
        identity.inchikey = p.get("InChIKey")
        ev = [Evidence(
            source="PubChem", kind="property",
            summary=f"{identity.formula}, MW {identity.molecular_weight}, CID {identity.pubchem_cid}",
            url=f"https://pubchem.ncbi.nlm.nih.gov/compound/{identity.pubchem_cid}")]
        return ev, SourceStatus(source="PubChem", ok=True, n_records=1)
    except Exception as e:
        return [], SourceStatus(source="PubChem", ok=False, error=str(e))


# ---------- PharmGKB / ClinPGx (drug-gene relationships) ----------
def pharmgkb(identity: DrugIdentity) -> tuple[list[Evidence], SourceStatus]:
    try:
        base = "https://api.clinpgx.org/v1/data"   # api.pharmgkb.org no longer resolves (PharmGKB -> ClinPGx)
        r = requests.get(f"{base}/chemical", params={"name": identity.query}, timeout=TIMEOUT, headers=HEADERS)
        r.raise_for_status()
        rows = r.json().get("data", [])
        if not rows:
            return [], SourceStatus(source="PharmGKB", ok=True, n_records=0, error="no chemical match")
        identity.pharmgkb_id = rows[0]["id"]
        identity.name = identity.name or rows[0].get("name")

        r = requests.get(f"{base}/clinicalAnnotation",
                         params={"relatedChemicals.accessionId": identity.pharmgkb_id, "view": "base"},
                         timeout=TIMEOUT, headers=HEADERS)
        r.raise_for_status()
        ev = []
        for a in r.json().get("data", []):
            loc = a.get("location", {}) or {}
            genes = [g.get("symbol") for g in loc.get("genes", []) if g.get("symbol")]
            level = (a.get("levelOfEvidence") or {}).get("term")
            types = ", ".join(a.get("types") or [])
            pops = None  # population detail needs a per-annotation fetch (TODO)
            ev.append(Evidence(
                source="PharmGKB", kind="pgx_annotation",
                gene=",".join(genes) or None,
                variant=loc.get("displayName"),
                summary=f"{types + ': ' if types else ''}{a.get('name', '')}".strip(),
                evidence_level=level, population=pops,
                url=f"https://www.clinpgx.org/clinicalAnnotation/{a.get('accessionId')}"))
        return ev, SourceStatus(source="PharmGKB", ok=True, n_records=len(ev))
    except Exception as e:
        return [], SourceStatus(source="PharmGKB", ok=False, error=str(e))


# ---------- Open Targets (targets, mechanisms, PGx) ----------
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

def _ot(query: str, variables: dict) -> dict:
    r = requests.post(OT_URL, json={"query": query, "variables": variables}, timeout=TIMEOUT, headers=HEADERS)
    r.raise_for_status()
    j = r.json()
    if j.get("errors"):
        raise RuntimeError(j["errors"][0].get("message", "GraphQL error"))
    return j["data"]


def opentargets(identity: DrugIdentity) -> tuple[list[Evidence], SourceStatus]:
    try:
        hits = _ot("""query($q:String!){search(queryString:$q, entityNames:["drug"]){hits{id name}}}""",
                   {"q": identity.query})["search"]["hits"]
        if not hits:
            return [], SourceStatus(source="OpenTargets", ok=True, n_records=0, error="no drug match")
        identity.chembl_id = hits[0]["id"]
        identity.name = identity.name or hits[0]["name"]

        d = _ot("""query($id:String!){drug(chemblId:$id){
                     id name drugType
                     mechanismsOfAction{rows{mechanismOfAction targets{approvedSymbol}}}
                   }}""", {"id": identity.chembl_id})["drug"]
        ev = []
        for row in (d.get("mechanismsOfAction") or {}).get("rows", []):
            genes = [t["approvedSymbol"] for t in row.get("targets") or []]
            ev.append(Evidence(
                source="OpenTargets", kind="mechanism", gene=",".join(genes) or None,
                summary=row.get("mechanismOfAction", ""),
                url=f"https://platform.opentargets.org/drug/{identity.chembl_id}"))

        # Pharmacogenomics field exists in newer schema versions; optional.
        try:
            pg = _ot("""query($id:String!){drug(chemblId:$id){pharmacogenomics{
                          variantRsId phenotypeText evidenceLevel target{approvedSymbol}}}}""",
                     {"id": identity.chembl_id})["drug"].get("pharmacogenomics") or []
            for x in pg:
                ev.append(Evidence(
                    source="OpenTargets", kind="pgx_annotation",
                    gene=(x.get("target") or {}).get("approvedSymbol"),
                    variant=x.get("variantRsId"), summary=x.get("phenotypeText") or "",
                    evidence_level=x.get("evidenceLevel"),
                    url=f"https://platform.opentargets.org/drug/{identity.chembl_id}"))
        except Exception:
            pass
        return ev, SourceStatus(source="OpenTargets", ok=True, n_records=len(ev))
    except Exception as e:
        return [], SourceStatus(source="OpenTargets", ok=False, error=str(e))


# ---------- PK-DB (pharmacokinetic studies) ----------
PKDB = "https://pk-db.com/api/v1"

def pkdb(identity: DrugIdentity) -> tuple[list[Evidence], SourceStatus]:
    """Lists PK studies for the substance. Parameter values (AUC, t1/2...) are not returned
    to anonymous users by /outputs/, so we report study metadata + subject counts only."""
    try:
        q = (identity.query or "").lower().strip()
        r = requests.get(f"{PKDB}/filter/", params={"studies__substance": q}, timeout=TIMEOUT, headers=HEADERS)
        r.raise_for_status()
        f = r.json()
        if not f.get("studies"):
            return [], SourceStatus(source="PK-DB", ok=True, n_records=0, error="substance not curated in PK-DB")
        r = requests.get(f"{PKDB}/studies/", params={"uuid": f["uuid"], "page_size": 50}, timeout=TIMEOUT, headers=HEADERS)
        r.raise_for_status()
        ev = []
        for st in r.json().get("data", {}).get("data", []):
            ref = st.get("reference") or {}
            ev.append(Evidence(
                source="PK-DB", kind="pk_study",
                summary=(f"{st.get('name')}: {ref.get('title') or ''} "
                         f"({st.get('individual_count', 0)} individuals, {st.get('group_count', 0)} groups, "
                         f"{st.get('output_count', 0)} PK outputs, {st.get('timecourse_count', 0)} time-courses)"),
                url=f"https://pk-db.com/data/{st.get('sid')}",
                variant=f"PMID {ref.get('pmid')}" if ref.get("pmid") else None))
        return ev, SourceStatus(source="PK-DB", ok=True, n_records=len(ev),
                                error=f"{f.get('individuals', 0)} individuals / {f.get('outputs', 0)} outputs indexed")
    except Exception as e:
        return [], SourceStatus(source="PK-DB", ok=False, error=str(e))


SOURCES = [pubchem, pharmgkb, opentargets, pkdb]
