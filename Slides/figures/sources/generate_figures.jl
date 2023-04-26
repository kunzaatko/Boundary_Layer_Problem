using CairoMakie
CairoMakie.activate!()

"""
    save_fig(name, fig)

Save a figure in both SVG and PDF formats.

# Parameters:
    - `name`: name of the file to save the figure as
    - `fig`: figure object to save

This function saves a figure object in both SVG and PDF formats. The SVG file is saved in the parent directory with the given name and the extension ".svg". The PDF file is saved in the same directory with the given name and the extension ".pdf". The Inkscape command is used to convert the SVG file to PDF format with text in LaTeX using PDF_TEX.
"""
function save_fig(name, fig; dpi=300)
    svg_path = "../" * name * ".svg"
    save(svg_path, fig)
    pdf_path = "../" * name * ".pdf"
    inkscape_cmd = Cmd(["inkscape", svg_path, "--export-area-page", "--export-dpi", string(dpi), "--export-type=pdf", "--export-latex", "--export-filename", pdf_path])
    run(inkscape_cmd)
end

for figure_path in filter(p -> splitext(p)[1] != "generate_figures" && splitext(p)[2] == ".jl", readdir())
    fig = include(figure_path)
    if fig isa AbstractVector
        for (name, f) in fig
            save_fig(name, f)
        end
    else
        save_fig(splitext(figure_path)[1], fig)
    end
end
