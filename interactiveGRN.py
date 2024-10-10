import marimo

__generated_with = "0.9.4"
app = marimo.App()


@app.cell
def __():
    import marimo
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy
    from scipy.integrate import solve_ivp
    import altair as alt
    import pandas as pd
    return alt, marimo, np, pd, plt, scipy, solve_ivp


@app.cell(hide_code=True)
def __(marimo):
    marimo.md("# Interactive exploration of biocircuits")
    return


@app.cell(hide_code=True)
def __(np, solve_ivp):
    # define a GRN system (simple example: A -> B)
    def grn_system(t, y, k1, k2, d1, d2, K):
        A, B = y
        dA_dt = k1 - d1 * A
        dB_dt = k2 * (A / (K + A)) - d2 * B
        return [dA_dt, dB_dt]

    # simulate a GRN system
    def simulate_grn(k1, k2, d1, d2, K, A0, B0,
                     t_start=0,
                     t_end=20,
                     num_points=500):
        initial_conditions = [A0, B0]
        t_eval = np.linspace(t_start, t_end, num_points)
        solution = solve_ivp(grn_system,
                             [t_start, t_end],
                             initial_conditions,
                             args=(k1, k2, d1, d2, K),
                             t_eval=t_eval)
        return solution.t, solution.y
    return grn_system, simulate_grn


@app.cell(hide_code=True)
def __(marimo):
    # Create sliders for parameters of the GRN (without 'description' argument)
    k1_slider = marimo.ui.slider(0, 5, value=1.0, step=0.1)
    k2_slider = marimo.ui.slider(0, 5, value=2.0, step=0.1)
    d1_slider = marimo.ui.slider(0, 2, value=0.5, step=0.05)
    d2_slider = marimo.ui.slider(0, 2, value=0.5, step=0.05)
    K_slider = marimo.ui.slider(0.1, 5.0, value=1.0, step=0.1)
    A0_slider = marimo.ui.slider(0, 5, value=0.0, step=0.1)
    B0_slider = marimo.ui.slider(0, 5, value=0.0, step=0.1)

    marimo.md(
        f"""
        Choose a value for $k_1$: {k1_slider} \n
        Choose a value for $k_2$: {k2_slider} \n
        Choose a value for $d_1$: {d1_slider} \n
        Choose a value for $d_2$: {d2_slider} \n
        Choose a value for $K$: {K_slider} \n
        Choose a value for $A_0$: {A0_slider} \n
        Choose a value for $B_0$: {B0_slider} \n
        """
    )
    return (
        A0_slider,
        B0_slider,
        K_slider,
        d1_slider,
        d2_slider,
        k1_slider,
        k2_slider,
    )


@app.cell(hide_code=True)
def __(
    A0_slider,
    B0_slider,
    K_slider,
    alt,
    d1_slider,
    d2_slider,
    k1_slider,
    k2_slider,
    np,
    pd,
    simulate_grn,
):
    # Extract slider values
    k1 = k1_slider.value
    k2 = k2_slider.value
    d1 = d1_slider.value
    d2 = d2_slider.value
    K = K_slider.value
    A0 = A0_slider.value
    B0 = B0_slider.value

    # Simulate the GRN system with the current parameters
    t, y = simulate_grn(k1, k2, d1, d2, K, A0, B0)

    # Extract the concentrations of A and B
    A = y[0]
    B = y[1]

    # Create a DataFrame to plot with Altair
    data = np.vstack([t, A, B]).T
    df = pd.DataFrame(data, columns=["Time", "Protein A", "Protein B"])

    # Create the Altair chart
    base = alt.Chart(df).encode(x='Time')

    line_A = base.mark_line(color='blue').encode(
        y='Protein A',
        tooltip=['Time', 'Protein A']
    )

    line_B = base.mark_line(color='green').encode(
        y='Protein B',
        tooltip=['Time', 'Protein B']
    )

    chart = alt.layer(line_A, line_B).properties(
        width=600,
        height=400,
        title='GRN dynamics over time'
    ).interactive()

    chart
    return (
        A,
        A0,
        B,
        B0,
        K,
        base,
        chart,
        d1,
        d2,
        data,
        df,
        k1,
        k2,
        line_A,
        line_B,
        t,
        y,
    )


if __name__ == "__main__":
    app.run()
