#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import streamlit as st

from peeling_core import (
    DEFAULTS,
    build_excel_bytes,
    make_inputs_from_display,
    result_summary,
    solve_auto,
)

st.set_page_config(
    page_title="Elastoplastic Peeling Analysis",
    page_icon="📐",
    layout="wide",
    initial_sidebar_state="collapsed",
)

st.markdown(
    """
    <style>
     /* 整个网页背景 */
    .stApp {
        background-color: #FFFFFF;
        color: #111827;
    }

    /* 主页面宽度 */
    .block-container {
        padding-top: 1.1rem;
        padding-bottom: 2rem;
        max-width: 1500px;
    }

    /* 主标题 */
    h1 {
        color: #17365D;
        font-weight: 800;
    }

    h2, h3, h4 {
        color: #244062;
    }

    /* 普通文字 */
    p, label, span {
        color: #111827;
    }

    /* 按钮 */
    div[data-testid="stButton"] > button,
    div[data-testid="stDownloadButton"] > button {
        min-height: 3rem;
        font-weight: 700;
        border-radius: 8px;
    }

    /* 输入框 */
    div[data-baseweb="input"] {
        background-color: #FFFFFF;
    }

    input {
        background-color: #FFFFFF !important;
        color: #111827 !important;
    }

    /* 结果卡片 */
    div[data-testid="stMetric"] {
        background: #FFFFFF;
        border: 2px solid #2F75B5;
        border-radius: 12px;
        padding: 14px 16px;
        box-shadow: 0 3px 10px rgba(0,0,0,0.08);
    }

    /* 结果名称 */
    div[data-testid="stMetricLabel"] {
        font-weight: 700;
        color: #44546A;
    }

    /* 结果数值 */
    div[data-testid="stMetricValue"] {
        font-size: 1.8rem;
        font-weight: 800;
        color: #005EB8;
    }

    /* Foundation model 区域 */
    .model-note {
        border: 1px solid #B4C7E7;
        background: #F4F8FC;
        color: #1F2937;
        border-radius: 10px;
        padding: .8rem 1rem;
        margin-bottom: .75rem;
    }

    /* Ready */
    .ready-box {
        background: #EFF8EF;
        border: 1px solid #70AD47;
        border-radius: 10px;
        padding: 1rem;
        color: #385723;
        font-family: monospace;
        font-weight: 700;
    }

    /* Detailed Output */
    textarea {
        background: #FFFFFF !important;
        color: #111827 !important;
    }

    </style>
    """,
    unsafe_allow_html=True,
)

# Intentionally no peel image in this web edition.
st.title("ELASTOPLASTIC PEELING ANALYSIS")
st.caption("Boundary match at x=0: y, y', y'' • No FPZ • V29 Codespaces edition")

for key, default in DEFAULTS.items():
    st.session_state.setdefault(key, default)

st.session_state.setdefault("result", None)
st.session_state.setdefault("excel_bytes", None)
st.session_state.setdefault("last_error", None)


def current_values():
    return {key: st.session_state[key] for key in DEFAULTS}


def load_example():
    for key, default in DEFAULTS.items():
        st.session_state[key] = default
    st.session_state.result = None
    st.session_state.excel_bytes = None
    st.session_state.last_error = None


def clear_output():
    st.session_state.result = None
    st.session_state.excel_bytes = None
    st.session_state.last_error = None


def calculate():
    try:
        film, adhesive, P, phi = make_inputs_from_display(current_values())
        result = solve_auto(film, adhesive, P, phi)
        summary = result_summary(result)

        st.session_state.result = result
        st.session_state.excel_bytes = build_excel_bytes(
            film, adhesive, P, phi, result, summary=summary
        )
        st.session_state.last_error = None

    except Exception as exc:
        st.session_state.result = None
        st.session_state.excel_bytes = None
        st.session_state.last_error = str(exc)


# Top command row, analogous to the V29 desktop layout.
_, c1, c2, c3, c4 = st.columns([5.0, 1.35, 1.35, 1.35, 1.35])

with c1:
    st.button("CALCULATE", type="primary", use_container_width=True, on_click=calculate)

with c2:
    st.button("LOAD EXAMPLE", use_container_width=True, on_click=load_example)

with c3:
    st.button("CLEAR OUTPUT", use_container_width=True, on_click=clear_output)

with c4:
    if st.session_state.excel_bytes:
        st.download_button(
            "EXPORT EXCEL",
            data=st.session_state.excel_bytes,
            file_name="peeling_results_v29_web.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            use_container_width=True,
        )
    else:
        st.button("EXPORT EXCEL", disabled=True, use_container_width=True)


left, right = st.columns([0.38, 0.62], gap="large")

with left:
    st.subheader("Foundation Model")
    st.markdown(
        """
        <div class="model-note">
        <b>Series foundation only:</b> k₁ (film) and k₂ (adhesive) are used in series.<br>
        No FPZ, σ<sub>f</sub>, δ<sub>f</sub>, or a<sub>f</sub> is used.<br>
        External transverse-force boundary: <b>V(0) = P sin(φ)</b>.
        </div>
        """,
        unsafe_allow_html=True,
    )

    with st.container(border=True):
        st.markdown("#### Film")
        st.number_input("Young's modulus E_f [GPa]", min_value=1e-6, step=0.05, format="%.6g", key="E_f_GPa")
        st.number_input("Poisson ratio nu_f [-]", min_value=-0.99, max_value=0.499, step=0.01, format="%.4f", key="nu_f")
        st.number_input("Yield stress sigma_y [MPa]", min_value=1e-6, step=1.0, format="%.6g", key="sigma_y_MPa")
        st.number_input("Hardening alpha [-]", min_value=0.0, max_value=0.999, step=0.01, format="%.4f", key="alpha")
        st.number_input("Width b_f [mm]", min_value=1e-6, step=0.5, format="%.6g", key="b_f_mm")
        st.number_input("Thickness h_f [um]", min_value=1e-6, step=1.0, format="%.6g", key="h_f_um")

    with st.container(border=True):
        st.markdown("#### Adhesive")
        st.number_input("Young's modulus E_a [GPa]", min_value=1e-6, step=0.1, format="%.6g", key="E_a_GPa")
        st.number_input("Poisson ratio nu_a [-]", min_value=-0.99, max_value=0.499, step=0.01, format="%.4f", key="nu_a")
        st.number_input("Width b_a [mm]", min_value=1e-6, step=0.5, format="%.6g", key="b_a_mm")
        st.number_input("Thickness h_a [um]", min_value=0.0, step=1.0, format="%.6g", key="h_a_um")

    with st.container(border=True):
        st.markdown("#### Peel Test")
        st.number_input("Peel force P [N]", min_value=0.0, step=0.05, format="%.6g", key="P_N")
        st.number_input("Peel angle phi [deg]", min_value=0.001, max_value=180.0, step=5.0, format="%.6g", key="phi_deg")

with right:
    st.subheader("Primary Result")

    if st.session_state.last_error:
        st.error(st.session_state.last_error)

    result = st.session_state.result

    if result is None:
        st.markdown('<div class="ready-box">Ready.</div>', unsafe_allow_html=True)
    else:
        r = result

        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric("theta_r", f"{r.peel.theta_deg:.6f}°")
        m2.metric("kappa_m", f"{r.peel.kappa_m:.6f}")
        m3.metric("K_m", f"{r.peel.K_m:.4e} 1/m")
        m4.metric("y(0)", f"{r.y0*1e6:.6f} μm")
        m5.metric("R", f"{r.peel.R:.6f} J/m²")

        st.markdown("#### Calculation Result")
        rows = [
            ("Model", r.name),
            ("theta_r [deg]", f"{r.peel.theta_deg:.9f}"),
            ("kappa_m", f"{r.peel.kappa_m:.12g}"),
            ("k_series [N/m²]", f"{r.k_series:.9e}"),
            ("beta [1/m]", f"{r.beta:.9e}"),
            ("1/beta [mm]", f"{r.length_scale*1e3:.9f}"),
            ("R [J/m²]", f"{r.peel.R:.9f}"),
            ("Force residual [N]", f"{r.force_residual:.3e}"),
            ("Curvature residual [1/m]", f"{r.residual:.3e}"),
        ]
        st.table({
            "Quantity": [row[0] for row in rows],
            "Value": [row[1] for row in rows],
        })

        st.markdown("#### Detailed Output")
        st.text_area(
            "Detailed calculation output",
            value=result_summary(r),
            height=570,
            disabled=True,
            label_visibility="collapsed",
        )

st.caption(
    "GitHub Codespaces web edition: no peel image and no local executable is required."
)
