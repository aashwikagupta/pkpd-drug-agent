"""Streamlit front end for the PK/PD drug-response agent.

Run:  streamlit run app.py
"""
from __future__ import annotations

import math
import os

import altair as alt
import pandas as pd
import streamlit as st

import agent
import ai
from pk_data import CLASS_INFO, DRUGS, GENE_PHENOTYPES, INTERACTIONS
from pkpd import HEPATIC, PatientProfile, flags, simulate

st.set_page_config(page_title="PK/PD Drug Response Agent", page_icon="💊", layout="wide")

# ------------------------------------------------------------------ styling
BLUE, GREY, BAND = "#2a78d6", "#8a8984", "#e9f1fb"
LEVEL = {  # status colors always paired with an icon + text label
    "critical": ("#d03b3b", "⛔", "Critical"),
    "warning": ("#c98500", "⚠️", "Caution"),
    "info": ("#2a78d6", "ℹ️", "Note"),
    "good": ("#1a8f5a", "✅", "OK"),
}
st.markdown("""
<style>
.block-container {padding-top: 1.6rem; max-width: 1350px;}
.hero h1 {font-size: 1.9rem; margin-bottom: 0;}
.hero p {color: #6b6a66; margin-top: .2rem;}
.pill {display:inline-block; padding:.12rem .55rem; border-radius:999px; font-size:.78rem;
       background:#f1f0ed; color:#3d3c39; margin-right:.3rem; border:1px solid #e3e2de;}
.flag {border:1px solid #e6e5e1; border-left-width:5px; border-radius:8px; padding:.55rem .8rem; margin-bottom:.5rem; background:#fff;}
.flag .t {font-weight:600; font-size:.93rem;}
.flag .d {color:#52514e; font-size:.86rem; margin-top:.15rem;}
.flag .tag {font-size:.72rem; color:#6b6a66; text-transform:uppercase; letter-spacing:.04em; margin-right:.4rem;}
.kpi {border:1px solid #e6e5e1; border-radius:10px; padding:.7rem .9rem; background:#fff; height:100%;}
.kpi .l {font-size:.78rem; color:#6b6a66;}
.kpi .v {font-size:1.55rem; font-weight:650; color:#0b0b0b; line-height:1.25;}
.kpi .s {font-size:.78rem; color:#6b6a66;}
.step {border:1px solid #e6e5e1; border-radius:8px; padding:.45rem .7rem; text-align:center; font-size:.82rem; background:#fafaf9;}
.step b {display:block; font-size:.9rem;}
@media (prefers-color-scheme: dark) {
  .flag, .kpi {background:#1f1f1e; border-color:#3a3a38;}
  .kpi .v {color:#f0efec;} .flag .d, .kpi .l, .kpi .s {color:#c3c2b7;}
  .pill, .step {background:#2a2a28; color:#e5e4df; border-color:#3a3a38;}
}
</style>
""", unsafe_allow_html=True)


def kpi(col, label, value, sub=""):
    col.markdown(f'<div class="kpi"><div class="l">{label}</div><div class="v">{value}</div>'
                 f'<div class="s">{sub}</div></div>', unsafe_allow_html=True)


def flag_card(f):
    color, icon, word = LEVEL[f.level]
    st.markdown(f'<div class="flag" style="border-left-color:{color}">'
                f'<div class="t">{icon} <span class="tag">{word} · {f.area}</span>{f.title}</div>'
                f'<div class="d">{f.detail}</div></div>', unsafe_allow_html=True)


def fmt(x, unit=""):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "–"
    s = f"{x:.3g}" if abs(x) < 1000 else f"{x:,.0f}"
    return f"{s} {unit}".strip()


# ------------------------------------------------------------------ sidebar: patient profile
with st.sidebar:
    st.markdown("### 💊 Drug & regimen")
    by_class: dict[str, list[str]] = {}
    for name, d in DRUGS.items():
        by_class.setdefault(d["cls"], []).append(name)
    options = [n for c in CLASS_INFO for n in by_class.get(c, [])]
    drug = st.selectbox("Drug", options, format_func=lambda n: f"{n.title()} ({DRUGS[n]['brand'].split(' ')[0]}) · {DRUGS[n]['cls']}")
    D = DRUGS[drug]
    daytime = "doses_per_day" in D
    c1, c2 = st.columns(2)
    dose = c1.number_input("Dose (mg)", min_value=0.05, max_value=3000.0, value=float(D["dose"]),
                           step=0.25 if D["dose"] < 5 else (5.0 if D["dose"] < 100 else 150.0), key=f"dose_{drug}")
    default_n = D.get("doses_per_day") or int(round(24 / D["interval"]))
    n_per_day = c2.selectbox("Doses / day", [1, 2, 3, 4], index=[1, 2, 3, 4].index(default_n), key=f"n_{drug}")
    interval = D["interval"] if daytime else 24 / n_per_day
    st.caption(f"Typical: {D['dose_range']}" + (" · daytime doses every 4 h" if daytime else f" · every {interval:g} h"))
    days_default = int(min(60, max(10, math.ceil(5 * D["t_half"] / 24) + 3)))
    days = st.slider("Simulate (days)", 2, 60, days_default, key=f"days_{drug}")
    adherence = st.slider("Adherence (% of doses taken)", 50, 100, 100, step=5) / 100

    st.markdown("### 🧍 Patient")
    c1, c2 = st.columns(2)
    age = c1.number_input("Age", 18, 100, 40)
    sex = c2.selectbox("Sex", ["Female", "Male"])
    weight = st.slider("Weight (kg)", 40, 150, 70)

    st.markdown("### 🩺 Organ function")
    egfr = st.slider("Kidney: eGFR (mL/min/1.73m²)", 10, 130, 100,
                     help="Matters most for renally cleared drugs (lithium).")
    hepatic = st.selectbox("Liver function", list(HEPATIC))

    st.markdown("### 🧬 Genotype → metabolizer phenotype")
    genes = list(dict.fromkeys(D["pgx_genes"] + [g for g, f in D["fm"].items() if g in GENE_PHENOTYPES and f >= 0.3]))
    phenotypes = {}
    if genes:
        for g in genes:
            choice = st.selectbox(g, ["Unknown / not tested"] + GENE_PHENOTYPES[g], key=f"pheno_{drug}_{g}")
            if choice != "Unknown / not tested":
                phenotypes[g] = choice
    else:
        st.caption("No CYP genotype modeled for this drug.")

    st.markdown("### 💊 Other medications & supplements")
    comeds = st.multiselect("Co-medications", [m for m in INTERACTIONS if m.lower() != drug],
                            help="Inhibitors/inducers change enzyme activity (phenoconversion). "
                                 "NSAIDs, ACE inhibitors, ARBs and thiazides change lithium clearance.")

    with st.expander("📝 Clinical context (optional, used by AI)"):
        cls = CLASS_INFO[D["cls"]]
        symptoms = st.text_area("Symptoms, diagnosis, treatment goals", placeholder="e.g. moderate depression with insomnia")
        prior = st.text_area("Previous medication responses", placeholder="e.g. nausea on fluoxetine")
        prefs = st.text_input("Lifestyle / routine / preferences", placeholder="e.g. night shifts, prefers once-daily")
        chars = {c: st.select_slider(c, ["Not assessed", "No concern", "Some concern", "Major concern"], key=f"ch_{drug}_{c}")
                 for c in cls["characteristics"]}

    st.markdown("### 🤖 AI reasoning (optional)")
    api_key = st.text_input("Anthropic API key", type="password",
                            help="Leave empty to use the ANTHROPIC_API_KEY environment variable. "
                                 "Without a key the agent still runs; only the AI narrative/chat is disabled.") or None
    st.caption("Research prototype. Not for clinical use.")

profile = PatientProfile(drug=drug, dose_mg=dose, interval_h=interval, doses_per_day=n_per_day if daytime else None,
                         days=days, age=age, sex=sex, weight_kg=weight, egfr=egfr, hepatic=hepatic,
                         phenotypes=phenotypes, comeds=comeds, adherence=adherence)
sim = simulate(profile)
flag_list = flags(profile, sim)
unit = sim.unit
m, mr = sim.metrics, sim.metrics_ref

# ------------------------------------------------------------------ chart helper
def conc_chart(s, rng, height=420):
    df = pd.DataFrame({"day": [t / 24 for t in s.t], "This patient": s.c_patient, "Typical adult": s.c_reference})
    step = max(1, len(df) // 1500)
    df = df.iloc[::step].reset_index(drop=True)
    long = df.melt("day", var_name="Series", value_name="conc")
    color = alt.Color("Series:N", scale=alt.Scale(domain=["This patient", "Typical adult"], range=[BLUE, GREY]),
                      legend=alt.Legend(orient="top", title=None))
    enc = dict(x=alt.X("day:Q", title="Time (days)"), y=alt.Y("conc:Q", title=f"Plasma concentration ({unit})"), color=color)
    # reference drawn first (dashed, behind), patient on top
    ref_line = alt.Chart(long[long.Series == "Typical adult"]).mark_line(strokeWidth=1.6, strokeDash=[5, 4]).encode(**enc)
    pat_line = alt.Chart(long[long.Series == "This patient"]).mark_line(strokeWidth=2.2).encode(**enc)
    lines = ref_line + pat_line
    layers = []
    if rng:
        band = pd.DataFrame({"lo": [rng[0]], "hi": [rng[1]], "label": [f"Reference range {rng[0]}–{rng[1]} {unit}"]})
        layers.append(alt.Chart(band).mark_rect(color=BAND, opacity=0.9).encode(y="lo:Q", y2="hi:Q"))
        layers.append(alt.Chart(band).mark_text(align="left", dx=4, dy=-6, color="#52514e", fontSize=11)
                      .encode(y="hi:Q", x=alt.value(0), text="label:N"))
    nearest = alt.selection_point(nearest=True, on="pointerover", fields=["day"], empty=False)
    sel = alt.Chart(df).mark_rule(color=GREY).encode(
        x="day:Q", opacity=alt.condition(nearest, alt.value(0.5), alt.value(0)),
        tooltip=[alt.Tooltip("day:Q", format=".2f", title="Day"),
                 alt.Tooltip("This patient:Q", format=".3g"), alt.Tooltip("Typical adult:Q", format=".3g")]
    ).add_params(nearest)
    pts = alt.Chart(long).mark_point(size=60, filled=True).encode(
        x="day:Q", y="conc:Q", color=color, opacity=alt.condition(nearest, alt.value(1), alt.value(0)))
    if s.missed_times:
        miss = pd.DataFrame({"day": [t / 24 for t in s.missed_times], "y": [0] * len(s.missed_times)})
        layers.append(alt.Chart(miss).mark_point(shape="triangle-up", size=70, color="#eb6834", filled=True)
                      .encode(x="day:Q", y="y:Q", tooltip=[alt.Tooltip("day:Q", format=".2f", title="Missed dose (day)")]))
    return alt.layer(*layers, lines, sel, pts).properties(height=height)


# ------------------------------------------------------------------ header
st.markdown(f'<div class="hero"><h1>PK/PD Drug Response Agent</h1>'
            f'<p>Patient-specific exposure and response for psychiatric drugs: pharmacokinetics, pharmacodynamics, '
            f'pharmacogenomics and live evidence in one report.</p></div>', unsafe_allow_html=True)
steps = st.columns(5)
for col, (a, b) in zip(steps, [("1 · Input", "drug + patient profile"), ("2 · Identify", "PubChem"),
                               ("3 · Retrieve", "ClinPGx · Open Targets · PK-DB"),
                               ("4 · Analyze", "PK sim · flags · conflicts"), ("5 · Output", "report · AI · CSV")]):
    col.markdown(f'<div class="step"><b>{a}</b>{b}</div>', unsafe_allow_html=True)
st.write("")

tabs = st.tabs(["📊 Overview", "📈 PK simulation", "🎯 PD & response", "🔎 Live evidence", "🤖 AI agent", "📄 Report"])
evidence_key = f"evidence_{drug}"
report = st.session_state.get(evidence_key)

# ------------------------------------------------------------------ Overview
with tabs[0]:
    st.markdown(f"#### {drug.title()} · {D['brand']}  "
                f"<span class='pill'>{D['cls']}</span><span class='pill'>Used for: {CLASS_INFO[D['cls']]['uses']}</span>",
                unsafe_allow_html=True)
    k = st.columns(5)
    ratio = sim.clearance.ratio
    kpi(k[0], "Exposure vs typical adult", f"{ratio:.2f}×", "AUC ratio from clearance")
    kpi(k[1], "Half-life", f"{sim.t_half_patient:.0f} h" if sim.t_half_patient < 100 else f"{sim.t_half_patient/24:.1f} d",
        f"typical {D['t_half']:g} h")
    kpi(k[2], "Time to steady state", f"{5 * sim.t_half_patient / 24:.1f} d", "≈ 5 half-lives")
    kpi(k[3], "Average level (last day)", fmt(m["cavg_ss"], unit), f"typical {fmt(mr['cavg_ss'], unit)}")
    if D.get("range"):
        kpi(k[4], "Time in reference range", f"{m.get('pct_in_range', 0):.0f}%", f"{D['range'][0]}–{D['range'][1]} {unit}, last 7 d")
    else:
        kpi(k[4], "Peak / trough (last day)", f"{fmt(m['cmax_ss'])} / {fmt(m['cmin_ss'])}", unit)
    st.write("")
    left, right = st.columns([1.15, 1])
    with left:
        st.markdown("##### What the agent found")
        for f in [f for f in flag_list if f.area != "Data"]:
            flag_card(f)
        with st.expander("Model assumptions & data caveats"):
            for f in [f for f in flag_list if f.area == "Data"]:
                flag_card(f)
    with right:
        st.markdown("##### Why exposure differs: clearance factors")
        fac = pd.DataFrame(sim.clearance.factors, columns=["Factor", "CL multiplier"])
        if fac.empty:
            st.info("No factor changes clearance from the typical 70 kg adult.")
        else:
            fac["Effect on exposure"] = 1 / fac["CL multiplier"]
            fac["Direction"] = fac["Effect on exposure"].map(lambda x: "Raises exposure" if x > 1 else "Lowers exposure")
            fac["log2"] = fac["Effect on exposure"].map(math.log2)
            span = max(1.0, fac["log2"].abs().max() * 1.1)
            bars = alt.Chart(fac).mark_bar(cornerRadiusEnd=4).encode(
                y=alt.Y("Factor:N", sort="-x", title=None, axis=alt.Axis(labelLimit=300)),
                x=alt.X("log2:Q", scale=alt.Scale(domain=[-span, span]),
                        axis=alt.Axis(labelExpr="format(pow(2, datum.value), '.2~f') + '×'", tickCount=5),
                        title="Effect on exposure (1× = no change)"),
                color=alt.Color("Direction:N", scale=alt.Scale(domain=["Raises exposure", "Lowers exposure"],
                                                               range=["#eb6834", BLUE]), legend=alt.Legend(orient="bottom", title=None)),
                tooltip=["Factor", alt.Tooltip("Effect on exposure:Q", format=".2f")]).properties(height=alt.Step(36))
            st.altair_chart(bars, width="stretch")
        st.markdown("##### Mini concentration curve")
        st.altair_chart(conc_chart(sim, D.get("range"), height=230), width="stretch")


with tabs[1]:
    st.markdown("#### Concentration over time")
    st.caption("One-compartment model, first-order absorption. Solid line: this patient. Dashed line: typical 70 kg adult, "
               "normal metabolizer, full adherence. " + ("Orange triangles = missed doses." if sim.missed_times else ""))
    st.altair_chart(conc_chart(sim, D.get("range")), width="stretch")

    c1, c2 = st.columns([1, 1])
    with c1:
        st.markdown("##### Steady-state metrics (last simulated day)")
        rows = [("Cmax", m["cmax_ss"], mr["cmax_ss"]), ("Cmin (trough)", m["cmin_ss"], mr["cmin_ss"]),
                ("Cavg", m["cavg_ss"], mr["cavg_ss"]), ("AUC 24 h", m["auc24_ss"], mr["auc24_ss"])]
        tbl = pd.DataFrame([(n, fmt(a, unit if "AUC" not in n else f"{unit}·h"), fmt(b, unit if "AUC" not in n else f"{unit}·h"),
                             f"{a / b:.2f}×" if b else "–") for n, a, b in rows],
                           columns=["Metric", "This patient", "Typical adult", "Ratio"])
        st.dataframe(tbl, hide_index=True, width="stretch")
    with c2:
        st.markdown("##### Model parameters")
        cl = sim.clearance
        params = pd.DataFrame([
            ("Clearance CL/F", f"{cl.cl_patient:.3g} L/h", f"{cl.cl_ref:.3g} L/h"),
            ("Volume V/F", f"{D['vd_f'] * weight:.0f} L", f"{D['vd_f'] * 70:.0f} L"),
            ("Half-life", f"{sim.t_half_patient:.1f} h", f"{D['t_half']:g} h"),
            ("Absorption ka (fit to Tmax)", f"{sim.ka:.2f} /h", f"Tmax {D['tmax']:g} h"),
            ("Accumulation ratio", f"{1 / (1 - math.exp(-math.log(2) / sim.t_half_patient * interval)):.2f}", ""),
        ], columns=["Parameter", "This patient", "Typical"])
        st.dataframe(params, hide_index=True, width="stretch")
    with st.expander("Enzyme pathways (fraction of clearance × activity in this patient)"):
        fm = D["fm"]
        if fm:
            st.dataframe(pd.DataFrame([(k, f"{v:.0%}", f"{cl.pathway_activity.get(k, 1):.2f}") for k, v in fm.items()],
                                      columns=["Pathway", "Fraction of CL (fm)", "Activity vs normal"]),
                         hide_index=True, width="stretch")
        st.caption(f"Renal fraction: {D['fe']:.0%}. Source: {D['source']}")
    curve = pd.DataFrame({"time_h": sim.t, f"patient_{unit}": sim.c_patient, f"typical_{unit}": sim.c_reference})
    st.download_button("⬇️ Download concentration-time CSV", curve.to_csv(index=False), f"{drug}_pk_curve.csv")

# ------------------------------------------------------------------ PD
with tabs[2]:
    cls = CLASS_INFO[D["cls"]]
    c1, c2 = st.columns([1.2, 1])
    with c1:
        st.markdown("#### Mechanism & targets")
        st.write(D["mechanism"])
        st.markdown(" ".join(f"<span class='pill'>🎯 {t}</span>" for t in D["targets"]), unsafe_allow_html=True)
        if report:
            mech = [e for e in report.evidence if e.kind == "mechanism"]
            if mech:
                st.caption("Open Targets (live): " + "; ".join(f"{e.summary} ({e.gene})" for e in mech))
        st.markdown("#### Time course: exposure vs effect")
        st.write(f"**PK steady state:** ~{5 * sim.t_half_patient / 24:.1f} days for this patient.  \n**PD onset:** {D['onset']}")
        st.markdown("#### Adverse effects to monitor")
        for a in D["adverse"]:
            st.markdown(f"- {a}")
        st.markdown("#### Genetics of response (PD)")
        st.write(D["pd_genes"])
    with c2:
        st.markdown(f"#### Patient characteristics to assess · {D['cls']}")
        st.caption("From the project notes: what a prescribing-support system would need to evaluate for this class.")
        for c, v in chars.items():
            icon = {"Not assessed": "◻️", "No concern": "✅", "Some concern": "⚠️", "Major concern": "⛔"}[v]
            st.markdown(f"{icon} **{c}**: {v}")
        st.markdown("#### PD-relevant flags")
        pd_flags = [f for f in flag_list if f.area in ("PD", "PGx")]
        for f in pd_flags or []:
            flag_card(f)
        if not pd_flags:
            st.caption("No PD-specific flags for this profile.")
        exp = sim.clearance.ratio
        st.markdown("#### Exposure → likely effect (qualitative)")
        if exp >= 1.4:
            st.write("Higher exposure: adverse-effect risk likely increased; benefit may not increase proportionally.")
        elif exp <= 0.7:
            st.write("Lower exposure: higher risk of non-response at this dose.")
        else:
            st.write("Exposure near typical: response expected to depend mostly on PD factors (diagnosis, targets, history).")

# ------------------------------------------------------------------ Live evidence
with tabs[3]:
    st.markdown("#### Live evidence from public databases")
    st.caption("PubChem (identity) · ClinPGx/PharmGKB (gene–drug annotations) · Open Targets (targets, PGx) · PK-DB (PK studies)")
    if st.button(f"🔎 Fetch evidence for {drug}", type="primary"):
        with st.spinner("Querying PubChem, ClinPGx, Open Targets, PK-DB…"):
            st.session_state[evidence_key] = agent.run(drug)
        st.rerun()
    if report:
        cols = st.columns(len(report.source_status))
        for col, s in zip(cols, report.source_status):
            kpi(col, "ClinPGx (PharmGKB)" if s.source == "PharmGKB" else s.source, f"{s.n_records}" if s.ok else "failed",
                ("✅ ok" if s.ok else "⛔ error") + (f" · {s.error[:60]}" if s.error else ""))
        i = report.identity
        st.markdown(f"<span class='pill'>PubChem CID {i.pubchem_cid}</span><span class='pill'>ChEMBL {i.chembl_id}</span>"
                    f"<span class='pill'>ClinPGx {i.pharmgkb_id}</span><span class='pill'>{i.formula} · MW {i.molecular_weight}</span>",
                    unsafe_allow_html=True)
        st.write("")
        df = agent.evidence_df(report)
        sub = st.tabs(["Gene–drug annotations", "PK studies (PK-DB)", "Conflicts & gaps", "All evidence"])
        with sub[0]:
            pg = df[df.kind == "pgx_annotation"]
            src = st.radio("Source", sorted(pg.source.unique()) or ["–"], horizontal=True)
            show = pg[pg.source == src][["gene", "variant", "evidence_level", "summary", "url"]]
            st.dataframe(show, hide_index=True, width="stretch",
                         column_config={"url": st.column_config.LinkColumn("link", display_text="open")})
            if not pg.empty:
                g = pg.dropna(subset=["gene"]).assign(gene=lambda x: x.gene.str.split(",")).explode("gene")
                counts = g.assign(gene=g.gene.str.strip()).groupby("gene").size() \
                    .sort_values(ascending=False).head(12).reset_index(name="annotations")
                st.markdown("**Most-annotated genes** (both sources)")
                st.bar_chart(counts.set_index("gene"), horizontal=True, color=BLUE,
                             height=60 + 26 * len(counts), x_label="Annotations", y_label="")
        with sub[1]:
            pk = df[df.kind == "pk_study"]
            if pk.empty:
                st.info(f"No PK-DB studies for {drug}. The simulator uses curated label/literature values (see the PK tab).")
            else:
                st.dataframe(pk[["variant", "summary", "url"]].rename(columns={"variant": "reference"}), hide_index=True,
                             width="stretch", column_config={"url": st.column_config.LinkColumn("study", display_text="open")})
                st.caption("PK-DB parameter values require a PK-DB account; study metadata is public.")
        with sub[2]:
            st.markdown("**Conflicts across sources** (gene reported by one source only)")
            st.write(", ".join(c.split(":")[0] for c in report.conflicts) or "none detected")
            st.markdown("**What remains uncertain / missing**")
            for g in report.gaps:
                st.markdown(f"- {g}")
        with sub[3]:
            st.dataframe(df, hide_index=True, width="stretch")
            st.download_button("⬇️ Evidence CSV", df.to_csv(index=False), f"{drug}_evidence.csv")
    else:
        st.info("Click the button to pull live data. Results are kept for this drug while the app is open.")

# ------------------------------------------------------------------ AI agent
full_report = agent.pkpd_markdown(profile, sim, flag_list, report)
context = full_report + "\n\n## Clinical context entered by user\n" + "\n".join(
    [f"- Symptoms/goals: {symptoms or 'not given'}", f"- Previous responses: {prior or 'not given'}",
     f"- Lifestyle/preferences: {prefs or 'not given'}"] + [f"- {c}: {v}" for c, v in chars.items()])

with tabs[4]:
    st.markdown("#### Ask the agent")
    if not ai.available(api_key):
        st.warning("No Anthropic API key found. The agent's analysis (simulation, flags, evidence, report) works without one; "
                   "add a key in the sidebar or set `ANTHROPIC_API_KEY` to enable the AI narrative and chat.")
        st.markdown("**Rule-based summary (no AI):**")
        for f in [f for f in flag_list if f.level in ("critical", "warning")][:6]:
            flag_card(f)
    else:
        hist_key = f"chat_{drug}"
        ctx_key = f"ctx_{drug}"
        if st.session_state.get(ctx_key) != context:          # profile changed -> new conversation
            st.session_state[ctx_key] = context
            st.session_state[hist_key] = []
        history = st.session_state[hist_key]
        c1, c2 = st.columns([1, 3])
        if c1.button("✨ Generate AI report", type="primary"):
            history.append({"role": "user", "content": "Write the structured drug response report for this patient, "
                                                        "answering the four research questions."})
        effort = c2.select_slider("Reasoning effort", ["low", "medium", "high"], value="medium", label_visibility="collapsed")
        for msg in history[:-1] if history and history[-1]["role"] == "user" else history:
            with st.chat_message(msg["role"]):
                st.markdown(msg["content"])
        q = st.chat_input("e.g. Why is exposure higher than typical? What if they stop omeprazole?")
        if q:
            history.append({"role": "user", "content": q})
        if history and history[-1]["role"] == "user":
            with st.chat_message("user"):
                st.markdown(history[-1]["content"])
            with st.chat_message("assistant"):
                answer = st.write_stream(ai.stream_answer(context, history, api_key, effort))
            history.append({"role": "assistant", "content": answer})
        st.caption(f"Model: {ai.MODEL}. Claude only explains the numbers computed by the agent; it is given the full report as context.")

# ------------------------------------------------------------------ Report
with tabs[5]:
    if not report:
        st.caption("Tip: fetch live evidence first to include ClinPGx annotations, PK-DB studies, conflicts and gaps.")
    st.markdown(full_report)
    st.download_button("⬇️ Download report (.md)", full_report, f"{drug}_drug_response_report.md")
