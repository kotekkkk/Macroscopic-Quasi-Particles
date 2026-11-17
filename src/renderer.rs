use plotters::prelude::*;

#[derive(Clone)]
pub struct ParticleFrame {
    pub positions: Vec<f64>,
    pub angles: Vec<f64>,
}

pub fn render(
    history: &[ParticleFrame],
    filename: &str,
    time: f64,
    fps: u32,
) -> Result<(), Box<dyn std::error::Error>> {

    if history.is_empty() || history[0].positions.is_empty() {
        return Ok(());
    }

    let n = history[0].positions.len();

    // Fixed dipole arrow length
    let arrow_length = 0.4;

    // Find maximum particle displacement
    let max_disp = history.iter()
        .flat_map(|frame| frame.positions.iter())
        .map(|x| x.abs())
        .fold(0.0_f64, f64::max);

    // Horizontal axis limits
    let margin = max_disp + arrow_length + 0.5;
    let x_min = -margin;
    let x_max = (n - 1) as f64 + margin;

    // Compensate for the chart's aspect ratio.
    // This keeps horizontal and vertical arrow components
    // visually proportional.
    let plot_width = 900.0;
    let plot_height = 500.0;

    let x_range = x_max - x_min;
    let y_limit = x_range * plot_height / (2.0 * plot_width);

    let root = BitMapBackend::gif(
        filename,
        (1000, 600),
        (1000.0/(fps as f64)) as u32,
    )?.into_drawing_area();

    // Determine how many frames to skip
    let stride = if time > 0.0 {
        (((1.0/(fps as f64)) / time).round() as usize).max(1)
    } else {
        1
    };

    for frame in history.iter().step_by(stride) {

        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .caption("Magnetoelastic Lattice", ("sans-serif", 30))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(50)
            .build_cartesian_2d(
                x_min..x_max,
                -y_limit..y_limit,
            )?;

        chart.configure_mesh()
            .x_desc("Lattice position")
            .y_desc("Dipole orientation")
            .draw()?;

        // Draw the equilibrium lattice
        chart.draw_series(LineSeries::new(
            (0..n).map(|i| (i as f64, 0.0)),
            &BLACK,
        ))?;

        // Draw particles and dipole arrows
        for (i, (&displacement, &theta)) in
            frame.positions.iter().zip(&frame.angles).enumerate()
        {
            // Actual particle position
            let x = i as f64 + displacement;

            // Dipole endpoint
            let dx = arrow_length * theta.cos();
            let dy = arrow_length * theta.sin();

            let tip = (x + dx, dy);

            // Arrowhead geometry
            let head_length = 0.12;
            let head_width = 0.06;

            let base_x = tip.0 - head_length * theta.cos();
            let base_y = tip.1 - head_length * theta.sin();

            let wing1 = (
                base_x - head_width * theta.sin(),
                base_y + head_width * theta.cos(),
            );

            let wing2 = (
                base_x + head_width * theta.sin(),
                base_y - head_width * theta.cos(),
            );

            // Draw particle
            chart.draw_series(std::iter::once(
                Circle::new((x, 0.0), 6, BLUE.filled())
            ))?;

            // Draw dipole shaft
            chart.draw_series(std::iter::once(
                PathElement::new(
                    vec![(x, 0.0), tip],
                    RED.stroke_width(3),
                )
            ))?;

            // Draw arrowhead
            chart.draw_series(std::iter::once(
                Polygon::new(
                    vec![tip, wing1, wing2],
                    RED.filled(),
                )
            ))?;
        }

        root.present()?;
    }

    Ok(())
}