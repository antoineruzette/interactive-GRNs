import marimo

__generated_with = "0.9.4"
app = marimo.App(
    app_title="interactive-repressilator",
    layout_file="layouts/repressilator.slides.json",
)


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
    marimo.md(r"""
    # Interactive Exploration of the Repressilator

    The **Repressilator** is a synthetic genetic network consisting of three proteins (A, B, C) that repress each other in a cyclic manner. The system is described by the following ordinary differential equations (ODEs):

    $$\frac{dA}{dt} = \frac{\alpha}{1 + C^n} - \beta A$$

    $$\frac{dB}{dt} = \frac{\alpha}{1 + A^n} - \beta B$$

    $$\frac{dC}{dt} = \frac{\alpha}{1 + B^n} - \beta C$$

    - \( \alpha \): Synthesis rate of each protein.
    - \( \beta \): Degradation rate of each protein.
    - \( n \): Hill coefficient, determining the cooperativity of repression.

    Elowitz, M., Leibler, S. A synthetic oscillatory network of transcriptional regulators. Nature 403, 335–338 (2000). [https://doi.org/10.1038/35002125](https://doi.org/10.1038/35002125)
    """)
    return


@app.cell(hide_code=True)
def __(np, solve_ivp):
    # Define the Repressilator system (A ->| B ->| C ->| A)
    def repressilator_system(t, y, alpha, beta, n):
        A, B, C = y
        dA_dt = alpha / (1 + C**n) - beta * A
        dB_dt = alpha / (1 + A**n) - beta * B
        dC_dt = alpha / (1 + B**n) - beta * C
        return [dA_dt, dB_dt, dC_dt]

    # Simulate the Repressilator system
    def simulate_repressilator(alpha, beta, n, A0, B0, C0, t_start=0, t_end=50, num_points=500):
        initial_conditions = [A0, B0, C0]
        t_eval = np.linspace(t_start, t_end, num_points)
        solution = solve_ivp(repressilator_system, [t_start, t_end], initial_conditions, args=(alpha, beta, n), t_eval=t_eval)
        return solution.t, solution.y
    return repressilator_system, simulate_repressilator


@app.cell(hide_code=True)
def __(marimo):
    # Create sliders for Repressilator parameters
    alpha_slider = marimo.ui.slider(0, 100, value=50, step=1)
    beta_slider = marimo.ui.slider(0, 10, value=1.0, step=0.1)
    n_slider = marimo.ui.slider(1, 5, value=2, step=1)
    A0_slider = marimo.ui.slider(0, 5, value=0.0, step=0.1)
    B0_slider = marimo.ui.slider(0, 5, value=0.0, step=0.1)
    C0_slider = marimo.ui.slider(0, 5, value=0.0, step=0.1)

    marimo.md(
        f"""
        Choose a value for $\\alpha$: {alpha_slider} \n
        Choose a value for $\\beta$: {beta_slider} \n
        Choose a value for $n$: {n_slider} \n
        Choose a value for $A_0$: {A0_slider} \n
        Choose a value for $B_0$: {B0_slider} \n
        Choose a value for $C_0$: {C0_slider} \n
        """
    )
    return (
        A0_slider,
        B0_slider,
        C0_slider,
        alpha_slider,
        beta_slider,
        n_slider,
    )


@app.cell(hide_code=True)
def __(
    A0_slider,
    B0_slider,
    C0_slider,
    alpha_slider,
    alt,
    beta_slider,
    n_slider,
    np,
    pd,
    simulate_repressilator,
):
    # Extract slider values
    alpha = alpha_slider.value
    beta = beta_slider.value
    n = n_slider.value
    A0 = A0_slider.value
    B0 = B0_slider.value
    C0 = C0_slider.value

    # Simulate the Repressilator system with the current parameters
    t, y = simulate_repressilator(alpha, beta, n, A0, B0, C0)

    # Extract the concentrations of A, B, and C
    A = y[0]
    B = y[1]
    C = y[2]

    # Create a DataFrame to plot with Altair
    data = np.vstack([t, A, B, C]).T
    df = pd.DataFrame(data, columns=["Time", "Protein A", "Protein B", "Protein C"])

    # Create the Altair chart for time vs concentration
    base_time = alt.Chart(df).encode(x='Time')

    line_A = base_time.mark_line(color='blue').encode(
        y='Protein A',
        tooltip=['Time', 'Protein A']
    )

    line_B = base_time.mark_line(color='green').encode(
        y='Protein B',
        tooltip=['Time', 'Protein B']
    )

    line_C = base_time.mark_line(color='red').encode(
        y='Protein C',
        tooltip=['Time', 'Protein C']
    )

    time_chart = alt.layer(line_A, line_B, line_C).properties(
        width=1000,
        height=400,
        title='Repressilator dynamics over time'
    ).interactive()

    # Create concentration vs concentration plots (A vs B, A vs C, B vs C)
    plot_A_vs_B = alt.Chart(df).mark_point().encode(
        x=alt.X('Protein A', scale=alt.Scale(zero=False)),
        y=alt.Y('Protein B', scale=alt.Scale(zero=False)),
        tooltip=['Protein A', 'Protein B']
    ).properties(title="A vs B")

    plot_A_vs_C = alt.Chart(df).mark_point().encode(
        x=alt.X('Protein A', scale=alt.Scale(zero=False)),
        y=alt.Y('Protein C', scale=alt.Scale(zero=False)),
        tooltip=['Protein A', 'Protein C']
    ).properties(title="A vs C")

    plot_B_vs_C = alt.Chart(df).mark_point().encode(
        x=alt.X('Protein B', scale=alt.Scale(zero=False)),
        y=alt.Y('Protein C', scale=alt.Scale(zero=False)),
        tooltip=['Protein B', 'Protein C']
    ).properties(title="B vs C")

    # Concatenate the concentration vs concentration plots horizontally
    conc_chart = alt.hconcat(plot_A_vs_B, plot_A_vs_C, plot_B_vs_C).resolve_scale(color='independent')

    # Combine the time-series chart and the concentration vs concentration plots
    combined_chart = alt.vconcat(time_chart, conc_chart)

    combined_chart
    return (
        A,
        A0,
        B,
        B0,
        C,
        C0,
        alpha,
        base_time,
        beta,
        combined_chart,
        conc_chart,
        data,
        df,
        line_A,
        line_B,
        line_C,
        n,
        plot_A_vs_B,
        plot_A_vs_C,
        plot_B_vs_C,
        t,
        time_chart,
        y,
    )


if __name__ == "__main__":
    app.run()
