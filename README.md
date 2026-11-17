# PK/PD Drug Response Agent (prototype)

Input: a psychiatric drug + a patient profile (age, weight, kidney/liver function, metabolizer genotype,
co-medications/supplements, adherence). Process: simulate exposure, flag unusual metabolism, pull live
evidence from public databases, compare sources. Output: a structured drug-response report (markdown + CSV),
with an optional AI narrative/chat.

Research prototype. Not for clinical use.

## Run
```
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
# optional: AI narrative + chat
export ANTHROPIC_API_KEY=...      # or paste a key in the sidebar
```

## Files
| File | Role |
|---|---|
| `app.py` | Streamlit front end: sidebar patient profile; tabs Overview, PK simulation, PD & response, Live evidence, AI agent, Report |
| `pk_data.py` | Curated PK/PD library for 8 drugs (one or more per class from the notes): sertraline, escitalopram, citalopram, fluoxetine, aripiprazole, lithium, alprazolam, methylphenidate. Includes enzyme pathways, CPIC/label guidance, reference ranges, interaction table |
| `pkpd.py` | Engine: covariate-adjusted clearance (genotype x inhibitors/inducers = phenoconversion, eGFR, liver, age, weight), one-compartment oral model, steady-state metrics, rule-based flags |
| `agent.py` | Coordinator: identity -> live sources -> normalize -> conflicts/gaps; `pkpd_markdown()` builds the full report |
| `sources.py` | Live API clients: PubChem, ClinPGx (PharmGKB), Open Targets, PK-DB |
| `ai.py` | Optional Claude step (claude-opus-5-5): explains the computed results, answers follow-ups. Never computes numbers itself |
| `models.py` | Shared pydantic schemas |

## How the agent answers the research questions
1. **ADME / how long it stays in the body**: simulated concentration-time curve, half-life, time to steady state, Cmax/Cmin/AUC vs a typical adult; PK-DB studies when available.
2. **Genetic variants**: metabolizer phenotype drives pathway activity; CPIC/label guidance shown per phenotype; live ClinPGx + Open Targets annotations.
3. **Demographic / ancestry differences**: noted in the report; live sources rarely return population detail (gap is flagged). Ancestry and ethnicity are treated as distinct.
4. **Other factors**: clearance breakdown chart shows how much age, eGFR, liver function, co-medications and adherence each change exposure.

## Known limits
- PK parameters are curated literature/label values (`pk_data.py`), approximate; verify against the cited source.
- One-compartment linear model; active metabolites (norfluoxetine, dehydro-aripiprazole, desmethylsertraline) not modeled.
- PD is qualitative (mechanism, onset, adverse effects, class characteristics); no quantitative receptor-occupancy model yet.
- PK-DB parameter values need an account (see ROADBLOCKS.md).
