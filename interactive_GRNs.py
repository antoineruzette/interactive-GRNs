import marimo

__generated_with = "0.9.4"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo
    import scipy
    import numpy
    import matplotlib
    from scipy.integrate import solve_ivp
    from marimo import App, Slider, Plot

    return App, Plot, Slider, matplotlib, mo, numpy, scipy, solve_ivp


@app.cell
def __(
    A0_slider,
    B0_slider,
    K_slider,
    d1_slider,
    d2_slider,
    k1_slider,
    k2_slider,
    plot,
):
    import numpy as np
    from scipy.integrate import solve_ivp

    def gene_regulatory_network(t, y, params):
        A, B = y
        k1, k2, d1, d2, K = params
        dA_dt = k1 - d1 * A
        dB_dt = k2 * (A / (K + A)) - d2 * B
        return [dA_dt, dB_dt]

    def simulate_gene_network(k1, k2, d1, d2, K, A0, B0, t_start=0, t_end=20, num_points=500):
        params = [k1, k2, d1, d2, K]
        initial_conditions = [A0, B0]
        t_eval = np.linspace(t_start, t_end, num_points)
        solution = solve_ivp(
            gene_regulatory_network,
            [t_start, t_end],
            initial_conditions,
            args=(params,),
            t_eval=t_eval,
            method='RK45'
        )
        return solution.t, solution.y

    def update_plot():
        # Get current slider values
        k1 = k1_slider.value
        k2 = k2_slider.value
        d1 = d1_slider.value
        d2 = d2_slider.value
        K  = K_slider.value
        A0 = A0_slider.value
        B0 = B0_slider.value

        # Run simulation
        t, y = simulate_gene_network(k1, k2, d1, d2, K, A0, B0)

        # Update the plot
        A = y[0]
        B = y[1]
        plot.update_data(t, A, B)

    return (
        gene_regulatory_network,
        np,
        simulate_gene_network,
        solve_ivp,
        update_plot,
    )


@app.cell
def __(App, Slider):
    app = App(title="Gene Regulatory Network Simulator")

    # Define sliders for each parameter
    k1_slider = Slider(min=0.0, max=5.0, value=1.0, step=0.1, description='k1 (Synthesis rate of A)')
    k2_slider = Slider(min=0.0, max=5.0, value=2.0, step=0.1, description='k2 (Max synthesis rate of B)')
    d1_slider = Slider(min=0.0, max=2.0, value=0.5, step=0.05, description='d1 (Degradation rate of A)')
    d2_slider = Slider(min=0.0, max=2.0, value=0.5, step=0.05, description='d2 (Degradation rate of B)')
    K_slider  = Slider(min=0.1, max=5.0, value=1.0, step=0.1, description='K (Activation constant)')

    # Initial conditions sliders
    A0_slider = Slider(min=0.0, max=5.0, value=0.0, step=0.1, description='Initial A')
    B0_slider = Slider(min=0.0, max=5.0, value=0.0, step=0.1, description='Initial B')

    return (
        A0_slider,
        B0_slider,
        K_slider,
        app,
        d1_slider,
        d2_slider,
        k1_slider,
        k2_slider,
    )


if __name__ == "__main__":
    app.run()
