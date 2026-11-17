"""PK/PD engine: patient covariates -> adjusted clearance -> concentration-time curve -> flags.

Model: one-compartment, first-order oral absorption, linear elimination, superposition of doses.
    C(t) = sum_i  D_i * ka / (V * (ka - ke)) * (exp(-ke*(t-t_i)) - exp(-ka*(t-t_i)))
Clearance is adjusted with a Rowland-Matin style model:
    CL = CL_ref * [ fe * renal_factor + sum_p fm_p * activity_p + f_other ] * hepatic * age
where activity_p combines genotype phenotype and co-medication inhibition/induction
(so a normal metabolizer on a strong inhibitor is "phenoconverted").

This is a teaching / research tool. Genotype alone cannot predict the complete PK profile.
"""
from __future__ import annotations

import math
import random
from typing import Optional

import numpy as np
from pydantic import BaseModel, Field

from pk_data import DRUGS, INTERACTIONS, PHENOTYPE_ACTIVITY


_trapz = getattr(np, "trapezoid", None) or np.trapz


class PatientProfile(BaseModel):
    drug: str
    dose_mg: float
    interval_h: float
    doses_per_day: Optional[int] = None      # None = around the clock (24/interval)
    days: int = 14
    age: int = 40
    sex: str = "Female"
    weight_kg: float = 70.0
    egfr: float = 100.0                      # mL/min/1.73m2
    hepatic: str = "Normal"                  # Normal | Mild (Child-Pugh A) | Moderate (B) | Severe (C)
    phenotypes: dict[str, str] = Field(default_factory=dict)
    comeds: list[str] = Field(default_factory=list)
    adherence: float = 1.0                   # fraction of doses taken
    seed: int = 7


class ClearanceBreakdown(BaseModel):
    cl_ref: float
    cl_patient: float
    ratio: float
    factors: list[tuple[str, float]]         # (label, multiplier on CL)
    pathway_activity: dict[str, float]
    phenoconversion: list[str]


class SimResult(BaseModel):
    t: list[float]
    c_patient: list[float]
    c_reference: list[float]
    dose_times: list[float]
    missed_times: list[float]
    unit: str
    t_half_ref: float
    t_half_patient: float
    ka: float
    metrics: dict[str, float]
    metrics_ref: dict[str, float]
    clearance: ClearanceBreakdown


# severity of hepatic impairment (0-1); multiplied by the drug's hepatic_sens = max fractional CL loss
HEPATIC = {"Normal": 0.0, "Mild (Child-Pugh A)": 0.6, "Moderate (Child-Pugh B)": 0.85, "Severe (Child-Pugh C)": 1.0}


def ka_from_tmax(ke: float, tmax: float) -> float:
    """Solve tmax = ln(ka/ke)/(ka-ke) for ka > ke by bisection."""
    f = lambda ka: math.log(ka / ke) / (ka - ke) - tmax
    lo, hi = ke * 1.0001, 50.0
    if f(hi) > 0:            # even very fast absorption can't reach tmax this small
        return hi
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def clearance(p: PatientProfile) -> ClearanceBreakdown:
    d = DRUGS[p.drug]
    v = d["vd_f"] * p.weight_kg
    cl_ref = math.log(2) / d["t_half"] * d["vd_f"] * 70.0     # reference adult = 70 kg
    # weight: allometric scaling of clearance
    w_fac = (p.weight_kg / 70.0) ** 0.75

    fe = d["fe"]
    fm = d["fm"]
    f_other = max(0.0, 1.0 - fe - sum(fm.values()))

    # pathway activity = genotype x co-medications
    activity, pheno = {}, []
    for path in fm:
        a = PHENOTYPE_ACTIVITY.get(p.phenotypes.get(path, "Normal metabolizer (NM)"), 1.0)
        for m in p.comeds:
            mult = INTERACTIONS.get(m, {}).get(path)
            if mult is not None:
                a *= mult
        activity[path] = a
        geno = p.phenotypes.get(path, "Normal metabolizer (NM)")
        if not geno.startswith("Poor") and a <= 0.25 and fm[path] >= 0.3:
            pheno.append(f"{path}: genotype {geno.split(' (')[0].lower()}, but co-medication drops activity to "
                         f"~{a:.0%} -> behaves like a poor metabolizer (phenoconversion).")

    renal = min(p.egfr / 100.0, 1.2)
    for m in p.comeds:
        r = INTERACTIONS.get(m, {}).get("renal_li")
        if r is not None and p.drug == "lithium":
            renal *= r

    hep = max(1.0 - d["hepatic_sens"] * HEPATIC.get(p.hepatic, 0.0), 0.1)
    age_fac = d["elderly"] if p.age >= 65 else 1.0

    nonrenal_metab = sum(fm[k] * activity[k] for k in fm) + f_other
    nonrenal_ref = sum(fm.values()) + f_other
    total = fe * renal + nonrenal_metab * hep * age_fac
    total_ref = fe + nonrenal_ref

    factors = [("Body weight (allometric)", w_fac)]
    if fe > 0:
        factors.append((f"Kidney function (eGFR {p.egfr:.0f}{', + interacting drugs' if renal != min(p.egfr/100, 1.2) else ''})",
                        (fe * renal + nonrenal_ref) / total_ref))
    for k in fm:
        if activity[k] != 1.0:
            factors.append((f"{k} activity x{activity[k]:.2f} (genotype + co-meds)",
                            (total_ref - fm[k] + fm[k] * activity[k]) / total_ref))
    if hep != 1.0:
        factors.append((f"Liver function ({p.hepatic})", (fe + (1 - fe) * hep) / total_ref))
    if age_fac != 1.0:
        factors.append(("Age >= 65", (fe + (1 - fe) * age_fac) / total_ref))

    factors = [(lbl, f) for lbl, f in factors if abs(f - 1.0) > 0.005]
    cl_p = cl_ref * w_fac * total / total_ref
    return ClearanceBreakdown(cl_ref=cl_ref, cl_patient=cl_p, ratio=cl_ref / cl_p,
                              factors=factors, pathway_activity=activity, phenoconversion=pheno)


def _dose_times(p: PatientProfile, days: int) -> list[float]:
    per_day = p.doses_per_day or max(1, int(round(24 / p.interval_h)))
    start = 8.0
    times = []
    for day in range(days):
        for k in range(per_day):
            t = day * 24 + start + k * p.interval_h
            if t < (day + 1) * 24 + start:
                times.append(t)
    return times


def _conc(t: np.ndarray, doses: list[float], amt: float, v: float, ka: float, ke: float) -> np.ndarray:
    c = np.zeros_like(t)
    coef = amt * ka / (v * (ka - ke))
    for td in doses:
        dt = t - td
        m = dt > 0
        c[m] += coef * (np.exp(-ke * dt[m]) - np.exp(-ka * dt[m]))
    return c


def _metrics(t, c, interval, days, d, per_day) -> dict[str, float]:
    last_day = (t >= (days - 1) * 24 + 8) & (t < days * 24 + 8)
    cl = c[last_day] if last_day.any() else c
    tl = t[last_day] if last_day.any() else t
    auc24 = float(_trapz(cl, tl)) if len(cl) > 1 else 0.0
    out = {"cmax_ss": float(cl.max()), "cmin_ss": float(cl.min()), "cavg_ss": auc24 / 24.0, "auc24_ss": auc24}
    rng = d.get("range")
    if rng:
        last_week = t >= max(0, (days - 7)) * 24
        inr = (c[last_week] >= rng[0]) & (c[last_week] <= rng[1])
        out["pct_in_range"] = float(inr.mean() * 100)
    return out


def simulate(p: PatientProfile) -> SimResult:
    d = DRUGS[p.drug]
    cb = clearance(p)
    v = d["vd_f"] * p.weight_kg
    v_ref = d["vd_f"] * 70.0
    ke_ref = math.log(2) / d["t_half"]
    ke = cb.cl_patient / v
    ka = ka_from_tmax(ke_ref, d["tmax"])
    if ka <= ke * 1.01:
        ka = ke * 1.5
    to_amt = d.get("mg_to_amount", 1.0)
    scale = 1000.0 if d["unit"] == "ng/mL" else 1.0     # mg/L -> ng/mL
    amt = p.dose_mg * to_amt

    days = p.days
    t = np.linspace(0, days * 24 + 8, int(min(6000, max(1500, days * 24 * 4))))
    planned = _dose_times(p, days)
    rnd = random.Random(p.seed)
    n_miss = int(round((1 - p.adherence) * len(planned)))
    missed = sorted(rnd.sample(planned, n_miss)) if n_miss else []
    taken = [x for x in planned if x not in set(missed)]

    c_pat = _conc(t, taken, amt, v, ka, ke) * scale
    c_ref = _conc(t, planned, amt, v_ref, ka, ke_ref) * scale
    per_day = p.doses_per_day or max(1, int(round(24 / p.interval_h)))
    return SimResult(
        t=t.tolist(), c_patient=c_pat.tolist(), c_reference=c_ref.tolist(),
        dose_times=taken, missed_times=missed, unit=d["unit"],
        t_half_ref=d["t_half"], t_half_patient=math.log(2) / ke, ka=ka,
        metrics=_metrics(t, c_pat, p.interval_h, days, d, per_day),
        metrics_ref=_metrics(t, c_ref, p.interval_h, days, d, per_day),
        clearance=cb,
    )


# ---------------------------------------------------------------- flags ("agent reasoning")
class Flag(BaseModel):
    level: str          # "critical" | "warning" | "info" | "good"
    area: str           # PK | PD | PGx | Interaction | Data
    title: str
    detail: str


def flags(p: PatientProfile, r: SimResult) -> list[Flag]:
    d = DRUGS[p.drug]
    out: list[Flag] = []
    u = r.unit
    fold = r.metrics["auc24_ss"] / r.metrics_ref["auc24_ss"] if r.metrics_ref["auc24_ss"] else 1.0
    fold_cl = r.clearance.ratio

    # exposure vs typical
    if fold_cl >= 2:
        out.append(Flag(level="critical", area="PK", title=f"Predicted exposure ~{fold_cl:.1f}x a typical adult",
                        detail="Clearance is markedly reduced by the combination of factors in the clearance breakdown. "
                               "Unusual metabolism: consider a lower dose, slower titration, or therapeutic drug monitoring."))
    elif fold_cl >= 1.4:
        out.append(Flag(level="warning", area="PK", title=f"Predicted exposure ~{fold_cl:.1f}x a typical adult",
                        detail="Moderately reduced clearance. Side effects may appear at standard doses."))
    elif fold_cl <= 0.6:
        out.append(Flag(level="warning", area="PK", title=f"Predicted exposure only ~{fold_cl:.2f}x a typical adult",
                        detail="Increased clearance (rapid metabolism or enzyme induction). Risk of non-response at standard doses."))
    else:
        out.append(Flag(level="good", area="PK", title=f"Exposure close to typical (~{fold_cl:.2f}x)",
                        detail="No major clearance-altering factor detected from the profile entered."))

    # therapeutic range
    rng = d.get("range")
    cav = r.metrics["cavg_ss"]
    if rng:
        if cav > rng[1]:
            out.append(Flag(level="warning", area="PK", title=f"Average steady-state level above reference range ({cav:.3g} > {rng[1]} {u})",
                            detail="AGNP reference ranges are population guides; confirm with a measured trough level."))
        elif cav < rng[0]:
            out.append(Flag(level="info", area="PK", title=f"Average steady-state level below reference range ({cav:.3g} < {rng[0]} {u})",
                            detail="May be sub-therapeutic for some patients; response, not level alone, should guide dosing."))
        if d.get("alert") and r.metrics["cmax_ss"] >= d["alert"]:
            out.append(Flag(level="critical", area="PK", title=f"Peak level reaches the laboratory alert level ({d['alert']} {u})",
                            detail="Toxicity risk. " + ("Lithium toxicity: tremor, confusion, ataxia. " if p.drug == "lithium" else "")))

    # steady state timing
    tss = 5 * r.t_half_patient
    if tss > p.days * 24:
        out.append(Flag(level="info", area="PK", title=f"Not yet at steady state after {p.days} days",
                        detail=f"~5 half-lives = {tss/24:.1f} days. Levels will keep rising; extend the simulation window to see the plateau."))
    if r.t_half_patient > 1.5 * r.t_half_ref:
        out.append(Flag(level="info", area="PK", title=f"Half-life prolonged: {r.t_half_patient:.0f} h vs {r.t_half_ref:.0f} h typical",
                        detail="Longer time to steady state and to washout; dose changes take longer to show their full effect."))

    # phenoconversion & PGx guidance
    for msg in r.clearance.phenoconversion:
        out.append(Flag(level="warning", area="Interaction", title="Phenoconversion", detail=msg))
    for gene, guide in d.get("cpic", {}).items():
        ph = p.phenotypes.get(gene)
        if ph and ph in guide:
            out.append(Flag(level="warning" if ph.startswith(("Poor", "Ultra")) else "info", area="PGx",
                            title=f"{gene} {ph}", detail=guide[ph]))

    # interactions
    for m in p.comeds:
        info = INTERACTIONS.get(m, {})
        affected = [k for k in d["fm"] if k in info and d["fm"][k] >= 0.1]
        if affected or (p.drug == "lithium" and "renal_li" in info):
            path = ", ".join(affected) if affected else "renal lithium clearance"
            out.append(Flag(level="warning", area="Interaction", title=f"{m}: {info.get('kind', '')}",
                            detail=f"Affects {path} for {p.drug}."))
        elif m in INTERACTIONS:
            out.append(Flag(level="info", area="Interaction", title=f"{m}: no major PK effect on {p.drug} modeled",
                            detail=f"{info.get('kind', '')}. Pharmacodynamic interactions (e.g. serotonergic, sedative) are not modeled."))
        if d["cls"].startswith("Antidepressant") and m in ("Fluoxetine", "Paroxetine", "Fluvoxamine", "St. John's wort"):
            out.append(Flag(level="critical", area="PD", title=f"Serotonin syndrome risk: {p.drug} + {m}",
                            detail="Two serotonergic agents. Combination is generally avoided."))

    # population / demographic factors from the notes
    if p.age >= 65:
        extra = {"citalopram": " FDA: max 20 mg/day over age 60 (QT).",
                 "alprazolam": " Benzodiazepines are on the Beers list: fall, fracture and cognitive risk.",
                 "sertraline": " Watch for hyponatremia.", "escitalopram": " Watch for hyponatremia; max 10 mg/day suggested."}
        out.append(Flag(level="warning" if p.drug in ("citalopram", "alprazolam") else "info", area="PD",
                        title="Older adult", detail="Reduced clearance and higher sensitivity." + extra.get(p.drug, "")))
    if p.drug == "citalopram" and p.dose_mg * (p.doses_per_day or 24 / p.interval_h) > 20 and \
            (p.age > 60 or p.phenotypes.get("CYP2C19", "").startswith("Poor")):
        out.append(Flag(level="critical", area="PD", title="Citalopram dose exceeds 20 mg/day cap",
                        detail="Label caps citalopram at 20 mg/day for age >60 or CYP2C19 PM because of QT prolongation."))
    if p.drug == "lithium":
        if p.egfr < 30:
            out.append(Flag(level="critical", area="PK", title="eGFR < 30: lithium generally avoided",
                            detail="Severe renal impairment: high toxicity risk, long half-life."))
        elif p.egfr < 60:
            out.append(Flag(level="warning", area="PK", title="eGFR < 60: reduce lithium dose and monitor levels closely",
                            detail="Lithium clearance falls roughly in proportion to GFR."))
    if p.adherence < 0.8:
        short = r.t_half_patient < 24
        out.append(Flag(level="warning" if short else "info", area="PK",
                        title=f"Adherence {p.adherence:.0%}: {len(r.missed_times)} missed doses in the simulation",
                        detail=("Short half-life: missed doses cause deep troughs (rebound / withdrawal / loss of effect)." if short
                                else "Long half-life buffers missed doses, but levels drift below the plateau.")))
    if p.hepatic != "Normal" and d["hepatic_sens"] > 0.3:
        out.append(Flag(level="warning", area="PK", title=f"Hepatic impairment ({p.hepatic})",
                        detail=f"{p.drug} is mainly cleared by the liver; label recommends lower or less frequent dosing."))

    # data provenance (the "what remains uncertain" part of the report)
    out.append(Flag(level="info", area="Data", title="Model assumptions",
                    detail="One-compartment model with curated literature parameters. Active metabolites, nonlinear kinetics, "
                           "food effects, and inter-individual variability beyond the entered factors are not modeled. "
                           "Genotype alone cannot predict the complete PK profile."))
    if d.get("notes"):
        out.append(Flag(level="info", area="Data", title="Drug-specific caveat", detail=d["notes"]))
    order = {"critical": 0, "warning": 1, "info": 2, "good": 3}
    return sorted(out, key=lambda f: order[f.level])
