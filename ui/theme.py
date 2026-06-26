"""Global visual theme for sc_tool.

Call apply_theme() at the top of every page (right after set_page_config).
Call next_page_button() at the bottom of a page when the current step is done.
"""

import streamlit as st

# ─── CSS ──────────────────────────────────────────────────────────────────────
_CSS = """
<style>

/* ── Fonts ──────────────────────────────────────────────────────────────── */
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap');
html, body, [class*="css"] {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
}

/* ── Background ─────────────────────────────────────────────────────────── */
[data-testid="stAppViewContainer"] {
    background: linear-gradient(145deg, #FDFCFF 0%, #F4F0FF 45%, #EDE8FF 100%);
    min-height: 100vh;
}
[data-testid="stMain"] {
    background: transparent;
}

/* ── Main content ───────────────────────────────────────────────────────── */
.main .block-container {
    padding-top: 2.5rem;
    padding-bottom: 5rem;
    max-width: 1100px;
}

/* ── Sidebar ────────────────────────────────────────────────────────────── */
[data-testid="stSidebar"] {
    background: rgba(255, 255, 255, 0.72) !important;
    backdrop-filter: blur(24px) saturate(180%);
    -webkit-backdrop-filter: blur(24px) saturate(180%);
    border-right: 1px solid rgba(124, 58, 237, 0.10);
    box-shadow: 4px 0 32px rgba(109, 40, 217, 0.06);
}
[data-testid="stSidebarNav"] {
    padding-top: 1rem;
}

/* ── Title ──────────────────────────────────────────────────────────────── */
h1 {
    background: linear-gradient(135deg, #7C3AED 0%, #4F46E5 100%);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
    font-weight: 800 !important;
    letter-spacing: -0.03em;
    line-height: 1.15 !important;
    margin-bottom: 0.25rem !important;
}

/* ── Subheaders ─────────────────────────────────────────────────────────── */
h2 {
    color: #3B1FAD !important;
    font-weight: 700 !important;
    letter-spacing: -0.02em;
    margin-top: 1.75rem !important;
}
h3 {
    color: #4C2BA0 !important;
    font-weight: 600 !important;
    letter-spacing: -0.01em;
}

/* ── Caption / small text ───────────────────────────────────────────────── */
[data-testid="stCaptionContainer"] p,
small {
    color: #7C6FB0 !important;
}

/* ── Metric cards ───────────────────────────────────────────────────────── */
[data-testid="metric-container"] {
    background: rgba(255, 255, 255, 0.65) !important;
    backdrop-filter: blur(14px) saturate(160%);
    -webkit-backdrop-filter: blur(14px) saturate(160%);
    border: 1px solid rgba(124, 58, 237, 0.10);
    border-radius: 18px !important;
    padding: 1.1rem 1.4rem !important;
    box-shadow: 0 2px 20px rgba(109, 40, 217, 0.07),
                inset 0 1px 0 rgba(255, 255, 255, 0.85);
    transition: transform 0.2s ease, box-shadow 0.2s ease;
}
[data-testid="metric-container"]:hover {
    transform: translateY(-2px);
    box-shadow: 0 8px 28px rgba(109, 40, 217, 0.13);
}
[data-testid="metric-container"] label {
    color: #7C6FB0 !important;
    font-size: 0.78rem !important;
    font-weight: 600 !important;
    text-transform: uppercase;
    letter-spacing: 0.06em;
}
[data-testid="metric-container"] [data-testid="stMetricValue"] {
    font-weight: 700 !important;
    color: #1E1B4B !important;
}

/* ── Primary buttons ────────────────────────────────────────────────────── */
[data-testid="baseButton-primary"] {
    background: linear-gradient(135deg, #7C3AED 0%, #6D28D9 100%) !important;
    color: #ffffff !important;
    border: none !important;
    border-radius: 12px !important;
    font-weight: 600 !important;
    font-size: 0.92rem !important;
    letter-spacing: 0.01em;
    box-shadow: 0 2px 14px rgba(109, 40, 217, 0.28) !important;
    transition: all 0.2s cubic-bezier(0.34, 1.56, 0.64, 1) !important;
    padding: 0.55rem 1.4rem !important;
}
[data-testid="baseButton-primary"]:hover {
    transform: translateY(-2px) scale(1.01) !important;
    box-shadow: 0 8px 24px rgba(109, 40, 217, 0.38) !important;
}
[data-testid="baseButton-primary"]:active {
    transform: translateY(0) scale(0.99) !important;
}

/* ── Secondary buttons ──────────────────────────────────────────────────── */
[data-testid="baseButton-secondary"] {
    background: rgba(255, 255, 255, 0.70) !important;
    color: #6D28D9 !important;
    border: 1px solid rgba(124, 58, 237, 0.20) !important;
    border-radius: 12px !important;
    font-weight: 500 !important;
    backdrop-filter: blur(8px);
    transition: all 0.2s ease !important;
}
[data-testid="baseButton-secondary"]:hover {
    background: rgba(237, 233, 254, 0.85) !important;
    border-color: rgba(124, 58, 237, 0.38) !important;
    transform: translateY(-1px) !important;
}

/* ── Alerts ─────────────────────────────────────────────────────────────── */
[data-testid="stAlert"] {
    border-radius: 14px !important;
    backdrop-filter: blur(10px);
    border-left-width: 3px !important;
}

/* ── Expander ───────────────────────────────────────────────────────────── */
[data-testid="stExpander"] {
    background: rgba(255, 255, 255, 0.58) !important;
    backdrop-filter: blur(12px);
    border: 1px solid rgba(124, 58, 237, 0.10) !important;
    border-radius: 14px !important;
    overflow: hidden;
    transition: box-shadow 0.2s ease;
}
[data-testid="stExpander"]:hover {
    box-shadow: 0 4px 20px rgba(109, 40, 217, 0.08);
}

/* ── Tabs ───────────────────────────────────────────────────────────────── */
[data-testid="stTabs"] [data-baseweb="tab-list"] {
    background: rgba(224, 215, 255, 0.45);
    border-radius: 14px;
    padding: 5px;
    gap: 4px;
    border: 1px solid rgba(124, 58, 237, 0.08);
}
[data-testid="stTabs"] [data-baseweb="tab"] {
    border-radius: 10px;
    font-weight: 500;
    color: #6B6B9E;
    transition: all 0.2s ease;
    padding: 0.45rem 1.1rem;
}
[data-testid="stTabs"] [aria-selected="true"][data-baseweb="tab"] {
    background: rgba(255, 255, 255, 0.92) !important;
    box-shadow: 0 2px 10px rgba(109, 40, 217, 0.12);
    color: #6D28D9 !important;
    font-weight: 600 !important;
}

/* ── Form ───────────────────────────────────────────────────────────────── */
[data-testid="stForm"] {
    background: rgba(255, 255, 255, 0.52);
    backdrop-filter: blur(14px);
    border: 1px solid rgba(124, 58, 237, 0.10);
    border-radius: 18px;
    padding: 1.5rem 1.5rem 1rem !important;
    box-shadow: 0 2px 20px rgba(109, 40, 217, 0.05);
}

/* ── Text inputs ────────────────────────────────────────────────────────── */
[data-baseweb="input"] > div {
    border-radius: 10px !important;
    background: rgba(255, 255, 255, 0.82) !important;
    border: 1px solid rgba(124, 58, 237, 0.14) !important;
    transition: border 0.15s ease, box-shadow 0.15s ease;
}
[data-baseweb="input"] > div:focus-within {
    border-color: #7C3AED !important;
    box-shadow: 0 0 0 3px rgba(124, 58, 237, 0.12) !important;
}

/* ── Selectbox ──────────────────────────────────────────────────────────── */
[data-baseweb="select"] > div {
    border-radius: 10px !important;
    background: rgba(255, 255, 255, 0.82) !important;
    border: 1px solid rgba(124, 58, 237, 0.14) !important;
}
[data-baseweb="select"] > div:focus-within {
    border-color: #7C3AED !important;
    box-shadow: 0 0 0 3px rgba(124, 58, 237, 0.12) !important;
}

/* ── Number input ───────────────────────────────────────────────────────── */
[data-testid="stNumberInput"] [data-baseweb="input"] > div {
    border-radius: 10px !important;
}

/* ── Plotly chart wrapper ───────────────────────────────────────────────── */
[data-testid="stPlotlyChart"] {
    background: rgba(255, 255, 255, 0.58);
    backdrop-filter: blur(12px);
    border-radius: 18px;
    border: 1px solid rgba(124, 58, 237, 0.08);
    padding: 0.5rem;
    box-shadow: 0 2px 20px rgba(109, 40, 217, 0.06);
}

/* ── DataFrame ──────────────────────────────────────────────────────────── */
[data-testid="stDataFrame"] {
    border-radius: 14px;
    overflow: hidden;
    border: 1px solid rgba(124, 58, 237, 0.08);
    box-shadow: 0 2px 16px rgba(109, 40, 217, 0.05);
}

/* ── Divider ────────────────────────────────────────────────────────────── */
hr {
    border: none !important;
    height: 1px !important;
    background: linear-gradient(
        90deg,
        transparent 0%,
        rgba(109, 40, 217, 0.18) 30%,
        rgba(109, 40, 217, 0.18) 70%,
        transparent 100%
    ) !important;
    margin: 2rem 0 !important;
}

/* ── Radio & checkbox ───────────────────────────────────────────────────── */
[data-testid="stRadio"] label,
[data-testid="stCheckbox"] label {
    font-weight: 500;
    color: #2D1B8A;
}

/* ── Navigation button section ──────────────────────────────────────────── */
.nav-section {
    margin-top: 2.5rem;
    padding-top: 1.5rem;
    border-top: 1px solid rgba(109, 40, 217, 0.12);
    display: flex;
    justify-content: flex-end;
}

</style>
"""


def apply_theme() -> None:
    """Inject the sc_tool glass UI theme. Call at the top of every page."""
    st.markdown(_CSS, unsafe_allow_html=True)


def next_page_button(label: str, page: str) -> None:
    """Render a right-aligned navigation button to the next page.

    Args:
        label: Display name of the destination page.
        page: File path relative to the app root (e.g. 'pages/02_QC.py').
    """
    st.markdown("---")
    _, col = st.columns([2, 1])
    with col:
        if st.button(f"Continue to {label}  →", type="primary", use_container_width=True):
            st.switch_page(page)
