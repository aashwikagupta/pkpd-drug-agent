# Roadblocks log

| # | Source/step | What broke | Guess at cause | Fix applied / proposed | Need from Sarah? |
|---|-------------|-----------|----------------|------------------------|------------------|
| 1 | PharmGKB | `api.pharmgkb.org` does not resolve (HTTP 000) | PharmGKB rebranded to ClinPGx | Switched to `api.clinpgx.org/v1/data` (same paths, works) | No |
| 2 | PK-DB | `/outputs/` returns 0 rows for every substance, even with a valid filter `uuid` | Parameter values hidden from anonymous users / some studies have closed licence | Use `/filter/?studies__substance=` + `/studies/?uuid=` for study metadata; simulator uses curated values | Do we have / can we get a PK-DB account? |
| 3 | PK-DB | No SSRIs, aripiprazole, lithium or methylphenidate curated; alprazolam (1 study), diazepam, midazolam, bupropion are | Coverage | Pick the "one selected drug" for PK-DB from what it covers (alprazolam?) | Which drug? |
| 4 | Open Targets | Works; `pharmacogenomics` field returns many annotations (e.g. 91 for sertraline) | - | - | No |
| 5 | Population data | Neither ClinPGx clinical annotations (base view) nor Open Targets PGx returns ancestry per finding | Needs per-annotation fetch / allele-frequency query | TODO: Open Targets variant allele frequencies (gnomAD) per rsID | - |

## Questions to bring to the meeting
- How many databases should the agent use for v1? (now 4: PubChem, ClinPGx, Open Targets, PK-DB)
- Which single drug for the PK-DB add-on? (PK-DB has alprazolam, not the SSRIs)
- ClinPGx access terms for programmatic use?
- Drug name matching: brand vs generic vs salt forms
- Should PD get a quantitative model (receptor occupancy / Emax) for one drug?
