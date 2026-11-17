"""Optional LLM step: Claude turns the agent's structured results into a narrative and answers follow-ups.

The deterministic pipeline (pkpd.py + agent.py) does all the numbers; Claude only explains them,
so every claim it makes can be traced back to the context block it was given.
Works only when an Anthropic API key is available; the app falls back to the rule-based report otherwise.
"""
from __future__ import annotations

import os
from typing import Iterator

import anthropic

MODEL = "claude-opus-5-5"

SYSTEM = """You are the reasoning layer of a research drug-response agent built by a student research team studying \
psychiatric drugs (pharmacokinetics, pharmacodynamics, pharmacogenomics).

You receive a CONTEXT block produced by deterministic code: a patient profile, simulated PK metrics, a clearance \
breakdown, rule-based flags, and evidence from PubChem, ClinPGx/PharmGKB, Open Targets and PK-DB.

Rules:
- Ground every number and claim in the CONTEXT. If something is not in it, say it is not available rather than guessing.
- Organize answers around the research questions: (1) ADME / how long the drug stays in the body, (2) genetic variants \
affecting metabolism or response, (3) evidence across demographic / ancestry groups (genetic ancestry and ethnicity are \
not interchangeable), (4) other factors (age, kidney/liver function, co-medications, adherence) explaining response or duration.
- Separate PK (exposure) from PD (effect, benefit, adverse-effect risk). Note that genotype alone cannot predict the full PK profile.
- End with a short "What remains uncertain" list.
- This is a research prototype, not clinical advice: do not tell anyone to start, stop or change a medication; \
phrase dose observations as what guidelines or the model suggest, for discussion with a clinician.
- Be concise; use markdown headings and bullets."""


def available(api_key: str | None = None) -> bool:
    return bool(api_key or os.environ.get("ANTHROPIC_API_KEY") or os.environ.get("ANTHROPIC_AUTH_TOKEN"))


def _client(api_key: str | None) -> anthropic.Anthropic:
    return anthropic.Anthropic(api_key=api_key) if api_key else anthropic.Anthropic()


def stream_answer(context: str, history: list[dict], api_key: str | None = None,
                  effort: str = "medium") -> Iterator[str]:
    """Yields text chunks. `history` is a list of {"role": "user"|"assistant", "content": str}, ending with a user turn."""
    client = _client(api_key)
    system = [
        {"type": "text", "text": SYSTEM},
        # context is stable across follow-up questions on the same profile -> cache it
        {"type": "text", "text": f"<CONTEXT>\n{context}\n</CONTEXT>", "cache_control": {"type": "ephemeral"}},
    ]
    try:
        with client.beta.messages.stream(
            model=MODEL,
            max_tokens=16000,
            system=system,
            messages=history,
            output_config={"effort": effort},
            # if the primary model declines, the API retries on a fallback model in the same call
            betas=["server-side-fallback-2026-07-01"],
            fallbacks="default",
        ) as stream:
            for text in stream.text_stream:
                yield text
            final = stream.get_final_message()
            if final.stop_reason == "refusal":
                yield "\n\n_The model declined to answer this request._"
            elif final.stop_reason == "max_tokens":
                yield "\n\n_(Answer truncated at the token limit.)_"
    except anthropic.AuthenticationError:
        yield "**Authentication failed.** Check the Anthropic API key in the sidebar."
    except anthropic.PermissionDeniedError as e:
        yield f"**Permission denied by the API:** {e.message}"
    except anthropic.RateLimitError:
        yield "**Rate limited.** Wait a moment and try again."
    except anthropic.BadRequestError as e:
        yield f"**Request rejected:** {e.message}"
    except anthropic.APIStatusError as e:
        yield f"**API error {e.status_code}:** {e.message}"
    except anthropic.APIConnectionError:
        yield "**Could not reach the Anthropic API.** Check your network connection."
