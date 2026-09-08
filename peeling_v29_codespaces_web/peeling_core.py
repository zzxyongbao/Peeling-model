#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

from dataclasses import dataclass
from io import BytesIO
from datetime import datetime
import math

@dataclass(frozen=True)
class Film:
    E: float
    nu: float
    sigma_y: float
    alpha: float
    b: float
    h: float

@dataclass(frozen=True)
class Adhesive:
    E: float
    nu: float
    b: float
    h: float

@dataclass
class PeelResult:
    theta: float
    theta_deg: float
    kappa_m: float
    K_m: float
    branch: str
    R: float
    external: float
    bending: float
    tensile: float
    M_root: float

@dataclass
class AutoResult:
    name: str
    peel: PeelResult
    beta: float
    length_scale: float
    K_foundation: float
    residual: float
    k1: float
    k2: float | None
    k_series: float
    k_used: float
    Eeq: float
    nueq: float
    Eprime_eq: float
    I_f: float
    D_att: float


def deg(x): return x * math.pi / 180.0
def rad2deg(x): return x * 180.0 / math.pi


def check_nu(nu, label):
    if not (-1.0 < nu < 0.5):
        raise ValueError(f"{label} Poisson ratio must satisfy -1 < nu < 0.5.")


def validate_film(f):
    if f.E <= 0 or f.sigma_y <= 0 or f.b <= 0 or f.h <= 0:
        raise ValueError("Film E, sigma_y, width and thickness must be positive.")
    check_nu(f.nu, "Film")
    if not (0 <= f.alpha < 1):
        raise ValueError("Hardening alpha must satisfy 0 <= alpha < 1.")


def validate_adhesive(a):
    if a.h < 0:
        raise ValueError("Adhesive thickness cannot be negative.")
    check_nu(a.nu, "Adhesive")
    if a.h > 0 and (a.E <= 0 or a.b <= 0):
        raise ValueError("Adhesive E and width must be positive when h_a > 0.")


# ---------- Latest foundation / attached-region model ----------

def eq12_k(E, nu, b, h):
    check_nu(nu, "Eq.12")
    if E <= 0 or b <= 0 or h <= 0:
        raise ValueError("Eq. (12) requires positive E, b and h.")
    return E * (1-nu) / ((1+nu)*(1-2*nu)) * (2*b/h)


def foundation_series(f, a):
    k1 = eq12_k(f.E, f.nu, f.b, f.h)
    if a.h == 0:
        return k1, None, k1
    k2 = eq12_k(a.E, a.nu, a.b, a.h)
    return k1, k2, k1*k2/(k1+k2)


def attached_properties(f, a):
    if a.h == 0:
        Eeq = f.E
        nueq = f.nu
    else:
        Eeq = (f.h+a.h) / (f.h/f.E + a.h/a.E)
        # Explicit engineering closure: thickness-weighted nu.
        nueq = (f.nu*f.h + a.nu*a.h) / (f.h+a.h)
    check_nu(nueq, "Equivalent")
    Eprime = Eeq / (1-nueq**2)
    I_f = f.b * f.h**3 / 12
    D_att = Eprime * I_f
    return Eeq, nueq, Eprime, I_f, D_att


# ---------- Chen peeled-arm elastoplastic mechanics ----------

def Me(f):
    return f.sigma_y*f.b*f.h**2 / (6*math.sqrt(1-f.nu+f.nu**2))


def Ke(f):
    return (
        f.sigma_y/(f.E*f.h)
        * 2*(1-f.nu**2)
        / math.sqrt(1-f.nu+f.nu**2)
    )


def eta(f, P):
    return 6*f.E*P/(f.sigma_y**2*f.b*f.h)


def reverse_threshold(alpha):
    return (2-alpha)/(1-alpha)


def m_loading(kappa, alpha):
    if kappa <= 1:
        return kappa
    return (1-alpha)/2*(3-1/kappa**2) + alpha*kappa


def eq5b_residual(k, et, angle_term, alpha):
    a = alpha
    return (
        a*((1-a)*k-(2-a))**2*k
        +(1-a)**2*(2-a)*k**2
        -2*(1-a)*(2-a)**2*k
        -et*k*angle_term
        +(2-a)**3
        +(2-a)*(k-1)*(a*k+2-a)*k
    )


def bisect(fun, lo, hi, tol=1e-11, max_iter=300):
    flo, fhi = fun(lo), fun(hi)
    if flo == 0: return lo
    if fhi == 0: return hi
    if flo*fhi > 0:
        raise ValueError("Root is not bracketed.")
    for _ in range(max_iter):
        mid = (lo+hi)/2
        fm = fun(mid)
        if abs(fm) < 1e-13 or abs(hi-lo) < tol*max(1, abs(mid)):
            return mid
        if flo*fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return (lo+hi)/2


def peak_kappa(f, P, phi, theta):
    angle_term = max(0.0, 1-math.cos(phi-theta))
    et = eta(f, P)
    trial = math.sqrt(max(0.0, et*angle_term))
    krev = reverse_threshold(f.alpha)

    if trial <= krev + 1e-12:
        return trial, ("elastic" if trial <= 1+1e-12 else "forward_plastic")

    fun = lambda k: eq5b_residual(k, et, angle_term, f.alpha)
    lo = krev
    upper = max(trial*1.5, krev+1)
    for _ in range(12):
        px, pf = lo, fun(lo)
        for i in range(1, 1000):
            x = lo + (upper-lo)*i/999
            fx = fun(x)
            if pf*fx <= 0:
                return bisect(fun, px, x), "reverse_plastic"
            px, pf = x, fx
        upper *= 2
    raise ValueError("Could not solve Chen Eq. (5b).")


def bending_integral(k, alpha):
    if k <= 1:
        return 0.0
    krev = reverse_threshold(alpha)
    if k <= krev:
        return (1-alpha)/2*(2/k + k**2 - 3)
    return 0.5 * (
        (1-alpha)*(
            alpha*(2-alpha)*k**2
            +3*(1-alpha)*(2-alpha)*k
            -3*(5-4*alpha+alpha**2)
        )
        +((2-alpha)**3 + 2*(1-alpha))/k
    )


def peel_result(f, P, phi, theta):
    km, branch = peak_kappa(f, P, phi, theta)
    Km = Ke(f) * km
    Mroot = Me(f) * m_loading(km, f.alpha)

    eps_t = P/(f.E*f.b*f.h)
    external = (P/f.b)*(1-math.cos(phi)+eps_t)
    bending = f.sigma_y**2*f.h/(3*f.E) * bending_integral(km, f.alpha)
    tensile = f.h*0.5*f.E*eps_t**2
    R = external - bending - tensile

    return PeelResult(
        theta, rad2deg(theta), km, Km, branch,
        R, external, bending, tensile, Mroot
    )


# ---------- Revised automatic theta_r ----------

def solve_auto(f, a, P, phi):
    """
    Solve theta_r and y(0) automatically.

    Attached decaying solution:
        y(x) = exp(beta*x)[A cos(beta*x) + B sin(beta*x)]

    At x=0:
        y(0)   = A
        y'(0)  = beta(A+B) = theta_r
        y''(0) = 2 beta^2 B = K_m(theta_r)

    Therefore:
        B  = K_m/(2 beta^2)
        y0 = theta_r/beta - K_m/(2 beta^2)

    The third derivative is:
        y'''(0) = 2 beta*K_m - 2 beta^2*theta_r

    External transverse-force equilibrium determines theta_r.

    Final v25 convention:
        V_peel = P sin(phi)

    This uses the same fixed global x-y coordinate system as the
    elastic-foundation beam equation.

    Chen Zhong's revised matching quantities remain y, y', and y''.
    The force equation is used as an external boundary condition,
    not as a fourth matching quantity.
    """
    validate_film(f)
    validate_adhesive(a)

    Eeq, nueq, Eprime, I_f, D_att = attached_properties(f, a)
    k1, k2, kseries = foundation_series(f, a)
    # v26 final foundation model: k1 and k2 are always in series.
    k_used = kseries
    name = "Series k1-k2 foundation"
    beta = (k_used/(4*D_att))**0.25

    def state(theta):
        p = peel_result(f, P, phi, theta)
        Km = p.K_m
        B = Km/(2*beta**2)
        y0 = theta/beta - B
        A = y0

        yp = beta*(A+B)
        ypp = 2*beta**2*B
        yppp = 2*beta**3*(B-A)

        V_att = D_att*yppp
        # Final v25 convention: global transverse-force component.
        V_peel = P*math.sin(phi)

        return p, A, B, y0, yp, ypp, yppp, V_att, V_peel, V_att-V_peel

    eps = max(1e-10, phi*1e-9)
    lo, hi = eps, phi-eps
    prev_t = lo
    prev_r = state(lo)[-1]
    brackets = []

    for i in range(1, 4001):
        th = lo + (hi-lo)*i/4000
        r = state(th)[-1]
        if math.isfinite(prev_r) and math.isfinite(r) and prev_r*r <= 0:
            brackets.append((prev_t, th))
        prev_t, prev_r = th, r

    if not brackets:
        raise ValueError(
            "No automatic root-angle solution was found from transverse-force "
            "equilibrium in 0 < theta_r < phi."
        )

    theta = bisect(lambda th: state(th)[-1], *brackets[0])
    p, A, B, y0, yp, ypp, yppp, V_att, V_peel, force_res = state(theta)

    result = AutoResult(
        name=name,
        peel=p,
        beta=beta,
        length_scale=1/beta,
        K_foundation=ypp,
        residual=p.K_m-ypp,
        k1=k1,
        k2=k2,
        k_series=kseries,
        k_used=k_used,
        Eeq=Eeq,
        nueq=nueq,
        Eprime_eq=Eprime,
        I_f=I_f,
        D_att=D_att,
    )

    result.A = A
    result.B = B
    result.y0 = y0
    result.slope_foundation = yp
    result.slope_residual = theta-yp
    result.yppp0 = yppp
    result.V_att = V_att
    result.V_peel = V_peel
    result.force_residual = force_res
    result.all_root_count = len(brackets)
    return result



DEFAULTS = {
    "E_f_GPa": 0.4,
    "nu_f": 0.40,
    "sigma_y_MPa": 100.0,
    "alpha": 0.05,
    "b_f_mm": 10.0,
    "h_f_um": 50.53,
    "E_a_GPa": 1.0,
    "nu_a": 0.44,
    "b_a_mm": 10.0,
    "h_a_um": 10.0,
    "P_N": 0.5,
    "phi_deg": 90.0,
}


def result_summary(result):
    r = result
    k2_text = "N/A (h_a=0)" if r.k2 is None else f"{r.k2:.9e} N/m²"

    return "\n".join([
        "REVISED CHEN ELASTIC-FOUNDATION MODEL — NO FPZ",
        "===============================================",
        "",
        "FOUNDATION",
        f"k1 film          = {r.k1:.9e} N/m²",
        f"k2 adhesive      = {k2_text}",
        f"k series         = {r.k_series:.9e} N/m²",
        "",
        "ATTACHED-REGION MODULUS",
        f"E_eq             = {r.Eeq/1e9:.9f} GPa",
        f"nu_eq            = {r.nueq:.9f}  [engineering estimate]",
        f"E'_eq            = {r.Eprime_eq/1e9:.9f} GPa",
        f"I_f              = {r.I_f:.9e} m^4",
        f"D_att            = {r.D_att:.9e} N·m²",
        "",
        "FULL MATCH AT x=0",
        f"y_att(0)=y_peel(0)  = {r.y0*1e6:.9f} um  [solved]",
        f"A = y(0)             = {r.A:.9e} m",
        f"B                     = {r.B:.9e} m",
        f"y'(0) foundation     = {r.slope_foundation:.9e} rad",
        f"theta_r peeled       = {r.peel.theta:.9e} rad = {r.peel.theta_deg:.9f} deg",
        f"slope residual        = {r.slope_residual:.9e} rad",
        f"y''(0) foundation    = {r.K_foundation:.9e} 1/m",
        f"K_m peeled           = {r.peel.K_m:.9e} 1/m",
        f"curvature residual   = {r.residual:.9e} 1/m",
        "",
        "EXTERNAL TRANSVERSE-FORCE BOUNDARY",
        "force boundary        = V = P sin(phi) [global]",
        f"third derivative y'''(0) = {r.yppp0:.9e} 1/m^2",
        f"V_att = D*y'''        = {r.V_att:.9e} N",
        f"V_peel                = {r.V_peel:.9e} N",
        f"force residual         = {r.force_residual:.9e} N",
        f"number of theta roots  = {r.all_root_count}",
        "",
        f"beta             = {r.beta:.9e} 1/m",
        f"1/beta           = {r.length_scale*1e3:.9f} mm",
        "",
        "PEELED ARM",
        f"kappa_m          = {r.peel.kappa_m:.12g}",
        f"branch           = {r.peel.branch}",
        f"M_root           = {r.peel.M_root:.9e} N·m",
        f"external work    = {r.peel.external:.9f} J/m²",
        f"bending work     = {r.peel.bending:.9f} J/m²",
        f"tensile work     = {r.peel.tensile:.9f} J/m²",
        f"R                = {r.peel.R:.9f} J/m²",
        "",
        "AUTOMATIC SOLVER",
        "1) y''(0) = K_m(theta_r)",
        "2) y'(0) = theta_r",
        "3) V_att(0) = V_peel determines theta_r",
        "4) y(0)=theta_r/beta-K_m/(2 beta^2) is then obtained",
    ])


def make_inputs_from_display(values):
    f = Film(
        float(values["E_f_GPa"]) * 1e9,
        float(values["nu_f"]),
        float(values["sigma_y_MPa"]) * 1e6,
        float(values["alpha"]),
        float(values["b_f_mm"]) * 1e-3,
        float(values["h_f_um"]) * 1e-6,
    )
    a = Adhesive(
        float(values["E_a_GPa"]) * 1e9,
        float(values["nu_a"]),
        float(values["b_a_mm"]) * 1e-3,
        float(values["h_a_um"]) * 1e-6,
    )
    P = float(values["P_N"])
    phi = deg(float(values["phi_deg"]))

    validate_film(f)
    validate_adhesive(a)
    if P < 0:
        raise ValueError("Peel force cannot be negative.")
    if not (0 < phi <= math.pi):
        raise ValueError("Peel angle must be in (0,180] deg.")
    return f, a, P, phi


def build_excel_bytes(f, a, P, phi, result, summary=None):
    import xlsxwriter

    if summary is None:
        summary = result_summary(result)

    stream = BytesIO()
    wb = xlsxwriter.Workbook(stream, {"in_memory": True})
    ws = wb.add_worksheet("Peeling Results")

    title_fmt = wb.add_format({
        "bold": True, "font_size": 16, "font_color": "#FFFFFF",
        "bg_color": "#173654", "align": "center", "valign": "vcenter"
    })
    section_fmt = wb.add_format({
        "bold": True, "font_color": "#FFFFFF",
        "bg_color": "#1F4E78", "border": 1
    })
    header_fmt = wb.add_format({
        "bold": True, "font_color": "#FFFFFF",
        "bg_color": "#4472C4", "border": 1, "align": "center"
    })
    label_fmt = wb.add_format({"bold": True, "bg_color": "#D9EAF7", "border": 1})
    value_fmt = wb.add_format({"border": 1, "num_format": "0.000000000"})
    sci_fmt = wb.add_format({"border": 1, "num_format": "0.000000E+00"})
    text_fmt = wb.add_format({"border": 1})
    wrap_fmt = wb.add_format({"text_wrap": True, "valign": "top", "border": 1})

    ws.merge_range("A1:H2", "Elastoplastic Peeling Analysis — V29 Web Results", title_fmt)
    ws.set_row(0, 22)
    ws.set_column("A:A", 31)
    ws.set_column("B:B", 20)
    ws.set_column("C:H", 18)

    row = 3
    ws.write(row, 0, "Run Information", section_fmt); row += 1
    info = [
        ("Timestamp", datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        ("Model", "V29 revised Chen elastic foundation — No FPZ"),
        ("Boundary matching", "y, y', y'' at x=0"),
        ("External force boundary", "V = P sin(phi)"),
        ("Root-angle mode", "Automatic only"),
    ]
    for label, val in info:
        ws.write(row, 0, label, label_fmt)
        ws.write(row, 1, val, text_fmt)
        row += 1

    row += 1
    ws.write(row, 0, "Input Parameters", section_fmt); row += 1
    inputs = [
        ("Film E_f [GPa]", f.E/1e9),
        ("Film nu_f [-]", f.nu),
        ("Film sigma_y [MPa]", f.sigma_y/1e6),
        ("Hardening alpha [-]", f.alpha),
        ("Film width b_f [mm]", f.b*1e3),
        ("Film thickness h_f [um]", f.h*1e6),
        ("Adhesive E_a [GPa]", a.E/1e9),
        ("Adhesive nu_a [-]", a.nu),
        ("Adhesive width b_a [mm]", a.b*1e3),
        ("Adhesive thickness h_a [um]", a.h*1e6),
        ("Peel force P [N]", P),
        ("Peel angle phi [deg]", rad2deg(phi)),
    ]
    for label, val in inputs:
        ws.write(row, 0, label, label_fmt)
        ws.write_number(row, 1, float(val), value_fmt)
        row += 1

    row += 1
    ws.write(row, 0, "Calculation Results", section_fmt); row += 1
    headers = [
        "Model", "theta_r [deg]", "kappa_m", "K_m [1/m]",
        "k [N/m^2]", "beta [1/m]", "1/beta [mm]", "R [J/m^2]"
    ]
    for col, h in enumerate(headers):
        ws.write(row, col, h, header_fmt)
    row += 1

    vals = [
        result.name, result.peel.theta_deg, result.peel.kappa_m,
        result.peel.K_m, result.k_used, result.beta,
        result.length_scale*1e3, result.peel.R
    ]
    ws.write(row, 0, vals[0], text_fmt)
    for c, val in enumerate(vals[1:], start=1):
        ws.write_number(row, c, float(val), sci_fmt if c in (3,4,5) else value_fmt)
    row += 1

    row += 1
    ws.write(row, 0, "Primary Series-Foundation Details", section_fmt); row += 1
    details = [
        ("k1 film [N/m^2]", result.k1),
        ("k2 adhesive [N/m^2]", result.k2),
        ("k_series [N/m^2]", result.k_series),
        ("E_eq [GPa]", result.Eeq/1e9),
        ("nu_eq [-]", result.nueq),
        ("E'_eq [GPa]", result.Eprime_eq/1e9),
        ("I_f [m^4]", result.I_f),
        ("D_att [N m^2]", result.D_att),
        ("Root displacement y(0) [um]", result.y0*1e6),
        ("Foundation constant A [m]", result.A),
        ("Foundation constant B [m]", result.B),
        ("Foundation slope y'(0) [rad]", result.slope_foundation),
        ("Slope residual [rad]", result.slope_residual),
        ("third derivative y'''(0) [1/m^2]", result.yppp0),
        ("V_att [N]", result.V_att),
        ("V_peel [N]", result.V_peel),
        ("Force residual [N]", result.force_residual),
        ("K_foundation [1/m]", result.K_foundation),
        ("Curvature residual [1/m]", result.residual),
        ("External work [J/m^2]", result.peel.external),
        ("Bending work [J/m^2]", result.peel.bending),
        ("Tensile work [J/m^2]", result.peel.tensile),
        ("M_root [N m]", result.peel.M_root),
    ]
    for label, val in details:
        ws.write(row, 0, label, label_fmt)
        if val is None:
            ws.write(row, 1, "N/A", text_fmt)
        else:
            num = float(val)
            use_sci = abs(num) >= 1e5 or (0 < abs(num) < 1e-4)
            ws.write_number(row, 1, num, sci_fmt if use_sci else value_fmt)
        row += 1

    row += 1
    ws.write(row, 0, "Detailed Output", section_fmt); row += 1
    ws.merge_range(row, 0, row+20, 7, summary, wrap_fmt)
    ws.set_row(row, 260)
    ws.freeze_panes(3, 0)

    wb.close()
    stream.seek(0)
    return stream.getvalue()


def self_test_web():
    f, a, P, phi = make_inputs_from_display(DEFAULTS)
    r = solve_auto(f, a, P, phi)

    assert abs(r.peel.theta_deg - 7.0275641405630696) < 1e-6
    assert abs(r.y0*1e6 - 1.7265221408994333) < 1e-6
    assert abs(r.peel.R - 50.061844448842265) < 1e-8
    assert abs(r.force_residual) < 1e-8
    assert abs(r.residual) < 1e-6

    xlsx = build_excel_bytes(f, a, P, phi, r)
    assert len(xlsx) > 5000 and xlsx[:2] == b"PK"

    return {
        "theta_r_deg": r.peel.theta_deg,
        "y0_um": r.y0*1e6,
        "R_J_m2": r.peel.R,
        "k_series_N_m2": r.k_series,
    }


if __name__ == "__main__":
    for key, value in self_test_web().items():
        print(f"{key} = {value}")
