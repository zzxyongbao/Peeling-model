# Elastoplastic Peeling Analysis — V29 Codespaces Web Edition

This package converts the V29 desktop program into a **browser-based Streamlit app** that runs in GitHub Codespaces.

## Main points

- The V29 calculation model is retained.
- The peel image is completely removed.
- No local `.exe` is required.
- Calculation, Load Example, Clear Output, and Excel export are retained.
- Excel files are generated in memory and downloaded through the browser.

## V29 model retained

- No FPZ.
- `k1` (film) and `k2` (adhesive) are used in series.
- Attached-region `E'` uses the film + adhesive series-modulus extension.
- Boundary matching at `x = 0`: `y`, `y'`, `y''`.
- External force boundary: `V(0) = P sin(phi)`.
- `theta_r` and `y(0)` are solved automatically.
- Chen peeled-arm elastoplastic calculation is retained for `kappa_m`, bending work and `R`.

## Default values

- Film `E_f = 0.4 GPa`
- Film `nu_f = 0.40`
- Film `sigma_y = 100 MPa`
- Hardening `alpha = 0.05`
- Film width `b_f = 10 mm`
- Film thickness `h_f = 50.53 um`
- Adhesive `E_a = 1.0 GPa`
- Adhesive `nu_a = 0.44`
- Adhesive width `b_a = 10 mm`
- Adhesive thickness `h_a = 10 um`
- Peel force `P = 0.5 N`
- Peel angle `phi = 90 deg`

## Run in GitHub Codespaces

1. Create a new GitHub repository.
2. Upload all files from this package, preserving the `.devcontainer` and `.streamlit` folders.
3. On GitHub select **Code → Codespaces → Create codespace on main**.
4. Wait while Codespaces builds the environment and installs `requirements.txt`.
5. The app starts automatically through `start.sh`.
6. Codespaces forwards port **8501** and normally opens the app in a browser.

If it does not open automatically:

- open the **PORTS** panel;
- find port **8501**;
- click **Open in Browser**.

## Manual start

```bash
bash start.sh
```

or:

```bash
streamlit run app.py
```

## Test the calculation core

```bash
python peeling_core.py
```

Expected V29 default result:

```text
theta_r ≈ 7.027564 deg
y(0)    ≈ 1.726522 um
R       ≈ 50.061844 J/m^2
```

## File structure

```text
app.py
peeling_core.py
requirements.txt
start.sh
README.md
.gitignore
.devcontainer/
  devcontainer.json
.streamlit/
  config.toml
```

There is intentionally **no peel image file**.
