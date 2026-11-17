"""Curated PK/PD library for the psychiatric drugs in the project notes.

Why curated: PubChem / PharmGKB / Open Targets carry no PK parameters, and PK-DB's
parameter values are not available to anonymous API users (see ROADBLOCKS.md). These
values are typical adult literature/label values used to drive the simulator. They are
APPROXIMATE and must be checked against the cited source before any real use.

Conventions
- vd_f      apparent volume of distribution (V/F), L/kg  -> oral dose is used as-is
- t_half    typical adult elimination half-life, h
- tmax      typical time to peak after an oral dose, h (absorption rate ka is fitted to it)
- fe        fraction of clearance that is renal (scaled by eGFR)
- fm        fraction of clearance through each enzyme pathway (Rowland-Matin model);
            whatever is left over is "other" and is not affected by genotype / inhibitors
- elderly   multiplier on non-renal clearance at age >= 65
- range     therapeutic reference range (AGNP consensus 2017 where available), same units
"""

# Enzyme activity relative to a normal metabolizer, applied to that pathway's fm.
PHENOTYPE_ACTIVITY = {
    "Poor metabolizer (PM)": 0.05,
    "Intermediate metabolizer (IM)": 0.6,
    "Normal metabolizer (NM)": 1.0,
    "Rapid metabolizer (RM)": 1.3,
    "Ultrarapid metabolizer (UM)": 1.6,
}
PHENOTYPES = list(PHENOTYPE_ACTIVITY)
# CYP2D6 has no "rapid" category in CPIC nomenclature.
GENE_PHENOTYPES = {
    "CYP2C19": PHENOTYPES,
    "CYP2D6": [p for p in PHENOTYPES if not p.startswith("Rapid")],
    "CYP2B6": [p for p in PHENOTYPES if not p.startswith("Ultrarapid")],
}

# Pathway activity multiplier for each interacting drug/supplement.
# Strong inhibitor ~0.2, moderate ~0.5, weak ~0.8; strong inducer ~3, moderate ~2.
# "renal_li" = reduces lithium renal clearance.
INTERACTIONS = {
    "Fluoxetine":            {"CYP2D6": 0.15, "kind": "strong CYP2D6 inhibitor (persists weeks after stopping)"},
    "Paroxetine":            {"CYP2D6": 0.15, "kind": "strong CYP2D6 inhibitor"},
    "Bupropion":             {"CYP2D6": 0.2,  "kind": "strong CYP2D6 inhibitor"},
    "Quinidine":             {"CYP2D6": 0.1,  "kind": "strong CYP2D6 inhibitor"},
    "Terbinafine":           {"CYP2D6": 0.2,  "kind": "strong CYP2D6 inhibitor"},
    "Fluvoxamine":           {"CYP2C19": 0.15, "CYP1A2": 0.15, "CYP3A4": 0.6, "kind": "strong CYP2C19/1A2 inhibitor"},
    "Fluconazole":           {"CYP2C19": 0.25, "CYP3A4": 0.5, "kind": "strong CYP2C19, moderate CYP3A4 inhibitor"},
    "Omeprazole":            {"CYP2C19": 0.5,  "kind": "moderate CYP2C19 inhibitor"},
    "Esomeprazole":          {"CYP2C19": 0.5,  "kind": "moderate CYP2C19 inhibitor"},
    "Ketoconazole":          {"CYP3A4": 0.1,  "kind": "strong CYP3A4 inhibitor"},
    "Itraconazole":          {"CYP3A4": 0.15, "kind": "strong CYP3A4 inhibitor"},
    "Clarithromycin":        {"CYP3A4": 0.2,  "kind": "strong CYP3A4 inhibitor"},
    "Ritonavir":             {"CYP3A4": 0.1,  "kind": "strong CYP3A4 inhibitor"},
    "Diltiazem":             {"CYP3A4": 0.5,  "kind": "moderate CYP3A4 inhibitor"},
    "Grapefruit juice":      {"CYP3A4": 0.6,  "kind": "dietary CYP3A4 inhibitor (supplement/diet)"},
    "Carbamazepine":         {"CYP3A4": 3.0, "CYP2C19": 2.0, "CYP2B6": 2.5, "kind": "strong CYP3A4 inducer"},
    "Rifampin":              {"CYP3A4": 4.0, "CYP2C19": 3.0, "CYP2B6": 3.0, "CYP2C9": 2.0, "kind": "strong pan-CYP inducer"},
    "Phenytoin":             {"CYP3A4": 2.5, "CYP2C19": 1.5, "kind": "strong CYP3A4 inducer"},
    "St. John's wort":       {"CYP3A4": 2.0, "CYP2C19": 1.5, "kind": "CYP3A4 inducer (herbal supplement)"},
    "Ibuprofen (NSAID)":     {"renal_li": 0.75, "kind": "NSAID: reduces lithium renal clearance"},
    "Naproxen (NSAID)":      {"renal_li": 0.75, "kind": "NSAID: reduces lithium renal clearance"},
    "Lisinopril (ACE-I)":    {"renal_li": 0.7,  "kind": "ACE inhibitor: reduces lithium renal clearance"},
    "Losartan (ARB)":        {"renal_li": 0.75, "kind": "ARB: reduces lithium renal clearance"},
    "Hydrochlorothiazide":   {"renal_li": 0.65, "kind": "thiazide diuretic: reduces lithium renal clearance"},
}

CLASS_INFO = {
    "Antidepressant (SSRI)": {
        "uses": "depression, anxiety disorders, PTSD, OCD",
        "characteristics": ["Anxiety level", "Sleep patterns", "Weight changes"],
    },
    "Anti-anxiety / sleep": {
        "uses": "anxiety and insomnia",
        "characteristics": ["Sleep patterns", "Sedation risk", "Fall risk"],
    },
    "Antipsychotic": {
        "uses": "psychosis, schizophrenia, bipolar symptoms",
        "characteristics": ["Symptom severity", "Metabolic health (weight, glucose, lipids)"],
    },
    "Mood stabilizer": {
        "uses": "bipolar disorder and mood-episode prevention",
        "characteristics": ["Mood-episode history", "Medical conditions (renal, thyroid)", "Drug interactions"],
    },
    "Stimulant": {
        "uses": "ADHD",
        "characteristics": ["Attention symptoms", "Sleep", "Appetite", "Blood pressure / heart rate"],
    },
}

DRUGS = {
    "sertraline": {
        "brand": "Zoloft", "cls": "Antidepressant (SSRI)",
        "dose": 50, "interval": 24, "dose_range": "25-200 mg/day",
        "vd_f": 40.0, "t_half": 26.0, "tmax": 6.0, "fe": 0.0,
        "fm": {"CYP2C19": 0.55, "CYP2B6": 0.2, "CYP3A4": 0.1, "CYP2D6": 0.05},
        "elderly": 0.7, "hepatic_sens": 0.6, "protein_binding": "98%",
        "unit": "ng/mL", "range": (10, 150), "alert": 300,
        "pgx_genes": ["CYP2C19", "CYP2B6"],
        "cpic": {
            "CYP2C19": {
                "Poor metabolizer (PM)": "CPIC: initiate a 50% lower starting dose and titrate to response, or choose an alternative not mainly metabolized by CYP2C19.",
                "Intermediate metabolizer (IM)": "CPIC: standard starting dose; consider slower titration and a lower maintenance dose.",
                "Ultrarapid metabolizer (UM)": "CPIC: standard starting dose; if no response at maintenance dose, consider an alternative.",
                "Rapid metabolizer (RM)": "CPIC: standard starting dose; if no response at maintenance dose, consider an alternative.",
            },
            "CYP2B6": {
                "Poor metabolizer (PM)": "CPIC: consider a lower starting dose / slower titration (CYP2B6 PM may raise sertraline exposure).",
            },
        },
        "mechanism": "Selective serotonin reuptake inhibitor: blocks the serotonin transporter (SERT, gene SLC6A4), raising synaptic serotonin.",
        "targets": ["SLC6A4"],
        "onset": "Some benefit in 1-2 weeks; full antidepressant effect typically 4-8 weeks (well after PK steady state at ~1 week).",
        "adverse": ["Nausea / GI upset", "Sexual dysfunction", "Insomnia or somnolence", "Activation / anxiety early on",
                    "Hyponatremia (older adults)", "Bleeding risk with NSAIDs / anticoagulants", "Serotonin syndrome with other serotonergic drugs"],
        "pd_genes": "Evidence for SLC6A4 (5-HTTLPR) and HTR2A variants predicting SSRI response is mixed and not used clinically.",
        "notes": "Active metabolite desmethylsertraline is ~10x less potent; not modeled.",
        "source": "Zoloft (sertraline) FDA label; CPIC SSRI guideline 2023 (Bousman et al.); Braten et al. 2020; AGNP TDM consensus 2017.",
    },
    "escitalopram": {
        "brand": "Lexapro", "cls": "Antidepressant (SSRI)",
        "dose": 10, "interval": 24, "dose_range": "10-20 mg/day",
        "vd_f": 20.0, "t_half": 30.0, "tmax": 5.0, "fe": 0.08,
        "fm": {"CYP2C19": 0.65, "CYP3A4": 0.25, "CYP2D6": 0.05},
        "elderly": 0.67, "hepatic_sens": 0.5, "protein_binding": "56%",
        "unit": "ng/mL", "range": (15, 80), "alert": 160,
        "pgx_genes": ["CYP2C19"],
        "cpic": {
            "CYP2C19": {
                "Poor metabolizer (PM)": "CPIC: consider a 50% lower starting dose and slower titration, or an alternative not mainly metabolized by CYP2C19.",
                "Intermediate metabolizer (IM)": "CPIC: standard starting dose; consider slower titration and lower maintenance dose.",
                "Rapid metabolizer (RM)": "CPIC: standard starting dose; if no response, consider titrating to a higher maintenance dose or an alternative.",
                "Ultrarapid metabolizer (UM)": "CPIC: consider an alternative not predominantly metabolized by CYP2C19 (low exposure expected).",
            },
        },
        "mechanism": "Selective serotonin reuptake inhibitor (S-enantiomer of citalopram); high-affinity SERT (SLC6A4) blockade.",
        "targets": ["SLC6A4"],
        "onset": "Full antidepressant / anxiolytic effect typically 4-8 weeks; steady state reached in ~1 week.",
        "adverse": ["Nausea", "Sexual dysfunction", "Insomnia / fatigue", "Dose-dependent QT prolongation", "Hyponatremia (older adults)"],
        "pd_genes": "SLC6A4 / HTR2A response associations are inconsistent; not actionable.",
        "notes": "",
        "source": "Lexapro FDA label; CPIC SSRI guideline 2023; Jukic et al. 2018 (Am J Psychiatry); AGNP 2017.",
    },
    "citalopram": {
        "brand": "Celexa", "cls": "Antidepressant (SSRI)",
        "dose": 20, "interval": 24, "dose_range": "20-40 mg/day (max 20 mg if >60 y or CYP2C19 PM)",
        "vd_f": 15.0, "t_half": 35.0, "tmax": 4.0, "fe": 0.12,
        "fm": {"CYP2C19": 0.6, "CYP3A4": 0.3, "CYP2D6": 0.05},
        "elderly": 0.75, "hepatic_sens": 0.55, "protein_binding": "80%",
        "unit": "ng/mL", "range": (50, 110), "alert": 220,
        "pgx_genes": ["CYP2C19"],
        "cpic": {
            "CYP2C19": {
                "Poor metabolizer (PM)": "CPIC: consider a 50% lower starting dose; FDA label caps the dose at 20 mg/day in CYP2C19 PMs (QT risk).",
                "Intermediate metabolizer (IM)": "CPIC: standard starting dose; consider slower titration and lower maintenance dose.",
                "Rapid metabolizer (RM)": "CPIC: standard starting dose; if no response, consider an alternative.",
                "Ultrarapid metabolizer (UM)": "CPIC: consider an alternative not predominantly metabolized by CYP2C19.",
            },
        },
        "mechanism": "Selective serotonin reuptake inhibitor (racemic); SERT (SLC6A4) blockade.",
        "targets": ["SLC6A4"],
        "onset": "Full effect typically 4-8 weeks.",
        "adverse": ["Dose-dependent QT prolongation (doses >40 mg not recommended)", "Nausea", "Sexual dysfunction", "Hyponatremia"],
        "pd_genes": "SLC6A4 / HTR2A response associations are inconsistent; not actionable.",
        "notes": "",
        "source": "Celexa FDA label (2011/2022 QT safety update); CPIC SSRI guideline 2023; AGNP 2017.",
    },
    "fluoxetine": {
        "brand": "Prozac", "cls": "Antidepressant (SSRI)",
        "dose": 20, "interval": 24, "dose_range": "20-80 mg/day",
        "vd_f": 30.0, "t_half": 96.0, "tmax": 7.0, "fe": 0.0,
        "fm": {"CYP2D6": 0.6, "CYP2C9": 0.2, "CYP3A4": 0.1},
        "elderly": 1.0, "hepatic_sens": 0.4, "protein_binding": "94.5%",
        "unit": "ng/mL", "range": None, "alert": None,
        "pgx_genes": ["CYP2D6"],
        "cpic": {
            "CYP2D6": {
                "Poor metabolizer (PM)": "CPIC 2023: no dosing recommendation for fluoxetine (insufficient evidence). Higher parent levels expected; active moiety (fluoxetine + norfluoxetine) changes less.",
            },
        },
        "mechanism": "Selective serotonin reuptake inhibitor; also a strong CYP2D6 inhibitor (perpetrator of interactions).",
        "targets": ["SLC6A4"],
        "onset": "Full effect 4-8 weeks; steady state takes ~4-5 weeks because of the very long half-life.",
        "adverse": ["Activation / insomnia", "Nausea", "Sexual dysfunction", "Long washout: interactions persist ~5 weeks after stopping"],
        "pd_genes": "SLC6A4 / HTR2A associations inconsistent.",
        "notes": "Parent drug only. Active metabolite norfluoxetine (t1/2 4-16 days) is not modeled; AGNP range (120-500 ng/mL) applies to fluoxetine + norfluoxetine, so no range is shown.",
        "source": "Prozac FDA label; CPIC SSRI guideline 2023; AGNP 2017.",
    },
    "aripiprazole": {
        "brand": "Abilify", "cls": "Antipsychotic",
        "dose": 10, "interval": 24, "dose_range": "10-30 mg/day",
        "vd_f": 5.6, "t_half": 75.0, "tmax": 4.0, "fe": 0.0,
        "fm": {"CYP2D6": 0.5, "CYP3A4": 0.45},
        "elderly": 0.8, "hepatic_sens": 0.9, "protein_binding": ">99%",
        "unit": "ng/mL", "range": (100, 350), "alert": 1000,
        "pgx_genes": ["CYP2D6"],
        "cpic": {
            "CYP2D6": {
                "Poor metabolizer (PM)": "FDA label: give half the usual dose (a quarter if also on a strong CYP3A4 inhibitor). DPWG: max 10 mg/day.",
            },
        },
        "mechanism": "Partial agonist at dopamine D2 (DRD2) and serotonin 5-HT1A (HTR1A) receptors; antagonist at 5-HT2A (HTR2A).",
        "targets": ["DRD2", "HTR1A", "HTR2A"],
        "onset": "Antipsychotic effect over 1-6 weeks; steady state ~2 weeks (~4 weeks in CYP2D6 PMs).",
        "adverse": ["Akathisia", "Weight gain / metabolic changes (lower than many antipsychotics)", "Nausea", "Impulse-control problems", "Tardive dyskinesia (long-term)"],
        "pd_genes": "DRD2 / HTR2A variants have exploratory associations with response and side effects; not actionable.",
        "notes": "Active metabolite dehydro-aripiprazole (~40% of parent at steady state) not modeled.",
        "source": "Abilify FDA label; DPWG aripiprazole guideline; AGNP 2017.",
    },
    "lithium": {
        "brand": "Lithobid / generic lithium carbonate", "cls": "Mood stabilizer",
        "dose": 300, "interval": 8, "dose_range": "900-1800 mg/day lithium carbonate (titrate to serum level)",
        "vd_f": 0.7, "t_half": 24.0, "tmax": 1.5, "fe": 0.95,
        "fm": {},
        "elderly": 1.0, "hepatic_sens": 0.0, "protein_binding": "0% (not protein bound)",
        "unit": "mmol/L", "mg_to_amount": 1 / 36.94,   # mg lithium carbonate -> mmol Li+
        "range": (0.6, 1.2), "alert": 1.5,
        "pgx_genes": [],
        "cpic": {},
        "mechanism": "Not fully understood: inhibits GSK-3 beta (GSK3B) and inositol monophosphatase (IMPA1); modulates second-messenger signalling.",
        "targets": ["GSK3B", "IMPA1"],
        "onset": "Antimanic effect 1-3 weeks; prophylaxis over months. Steady state ~5 days (longer with low eGFR).",
        "adverse": ["Narrow therapeutic index: toxicity >1.5 mmol/L (tremor, confusion, ataxia)", "Polyuria / nephrogenic diabetes insipidus",
                    "Hypothyroidism", "Long-term renal impairment", "Weight gain"],
        "pd_genes": "No actionable PGx; exploratory GWAS loci for response (e.g. ConLiGen) not clinical.",
        "notes": "Eliminated almost entirely by the kidney, so eGFR, dehydration, NSAIDs, ACE inhibitors/ARBs and thiazides dominate exposure. Concentrations in mmol/L Li+.",
        "source": "Lithium carbonate FDA label; AGNP 2017; ISBD/IGSLi lithium level consensus.",
    },
    "alprazolam": {
        "brand": "Xanax", "cls": "Anti-anxiety / sleep",
        "dose": 0.5, "interval": 8, "dose_range": "0.25-4 mg/day (divided)",
        "vd_f": 0.9, "t_half": 11.2, "tmax": 1.5, "fe": 0.2,
        "fm": {"CYP3A4": 0.75},
        "elderly": 0.7, "hepatic_sens": 0.55, "protein_binding": "80%",
        "unit": "ng/mL", "range": None, "alert": None,
        "pgx_genes": [],
        "cpic": {},
        "mechanism": "Positive allosteric modulator of GABA-A receptors (benzodiazepine site, e.g. GABRA1/GABRG2).",
        "targets": ["GABRA1", "GABRG2"],
        "onset": "Anxiolytic effect within ~1 hour of a dose.",
        "adverse": ["Sedation / psychomotor impairment", "Falls (older adults)", "Dependence, withdrawal and rebound anxiety (short half-life)",
                    "Respiratory depression with opioids or alcohol"],
        "pd_genes": "No actionable PGx.",
        "notes": "Short half-life: missed doses cause inter-dose withdrawal / rebound anxiety.",
        "source": "Xanax FDA label (ketoconazole contraindicated; elderly t1/2 16.3 h; alcoholic liver disease t1/2 19.7 h).",
    },
    "methylphenidate": {
        "brand": "Ritalin (IR)", "cls": "Stimulant",
        "dose": 10, "interval": 4, "doses_per_day": 3, "dose_range": "10-60 mg/day IR (2-3 doses)",
        "vd_f": 10.0, "t_half": 2.5, "tmax": 2.0, "fe": 0.0,
        "fm": {"CES1": 0.8},
        "elderly": 1.0, "hepatic_sens": 0.3, "protein_binding": "10-33%",
        "unit": "ng/mL", "range": None, "alert": None,
        "pgx_genes": [],
        "cpic": {},
        "mechanism": "Blocks dopamine (SLC6A3/DAT) and norepinephrine (SLC6A2/NET) reuptake.",
        "targets": ["SLC6A3", "SLC6A2"],
        "onset": "Effect within 30-60 min; lasts ~3-4 h (IR). No accumulation with standard dosing.",
        "adverse": ["Appetite suppression / weight loss", "Insomnia", "Increased heart rate and blood pressure", "Anxiety / irritability"],
        "pd_genes": "CES1 loss-of-function variants (e.g. G143E) raise exposure; evidence limited, not in CPIC guidance.",
        "notes": "Modeled as IR tablet; dosing interval set to 4 h for a daytime schedule. Extended-release forms differ.",
        "source": "Ritalin FDA label; PharmGKB CES1-methylphenidate annotations.",
    },
}

# Interval choices shown in the UI (hours)
INTERVALS = {"Every 4 h": 4, "Every 8 h (TID)": 8, "Every 12 h (BID)": 12, "Every 24 h (daily)": 24}
