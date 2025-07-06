import pandas as pd
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
from plotly.utils import PlotlyJSONEncoder
import json
import re
import os
import numpy as np
import argparse
import pysam
from collections import defaultdict


class GenomeCoverage:

    def __init__(self, mapping_file, sample_name, bed_file, bam_file, species, bin_size, output_dir):
        self.mapping_file = mapping_file
        self.sample_name = sample_name
        self.bed_file = bed_file
        self.bam_file = bam_file
        self.species = species
        self.bin_size = bin_size
        self.output_html = os.path.join(
            output_dir,
            f"coverage_report_{species.replace(' ', '_')}_{bin_size // 1000}kb.html"
        )
        os.makedirs(output_dir, exist_ok=True)
        self.df = self.read_and_merge()
        self.read_counts_per_species = self.count_mapped_reads_by_species()
        self.regions_sorted = sorted(self.df["region"].unique(), key=self.region_sort_key)

    def count_mapped_reads_by_species(self):
        print(f"Reading BAM file {self.bam_file} to count mapped reads...")
        species_per_chrom = {}

        # Lade die Mapping-Datei (chrom → species)
        mapping_df = pd.read_csv(self.mapping_file, sep="\t", names=["chrom", "species", "region"])
        for _, row in mapping_df.iterrows():
            species_per_chrom[row["chrom"]] = row["species"]

        # Zähle gemappte Reads pro Species
        bam = pysam.AlignmentFile(self.bam_file, "rb")
        counts = defaultdict(int)

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue  # nur primäre, gemappte Reads zählen
            ref = bam.get_reference_name(read.reference_id)
            species = species_per_chrom.get(ref)
            if species:
                counts[species] += 1

        bam.close()
        return dict(counts)

    def read_and_merge(self):
            print(f"Reading chromosome mapping file {self.mapping_file}...")
            mapping_df = pd.read_csv(self.mapping_file, sep="\t", names=["chrom", "species", "region"])
            print(f"Reading coverage BED file {self.bed_file}...")
            bed_df = pd.read_csv(self.bed_file, sep="\t", header=None, names=["chrom", "start", "end", "coverage"])
            bed_df["length"] = bed_df["end"] - bed_df["start"]
            merged_df = bed_df.merge(mapping_df, on="chrom")
            return merged_df[merged_df["species"] == self.species].copy()

    @staticmethod
    def region_sort_key(region):
        region = region.lower()
        if "chromosome" in region:
            match = re.search(r"chromosome\s+(\d+|x|y|z|w)", region)
            if match:
                val = match.group(1).lower()
                if val.isdigit():
                    return (0, int(val))
                return (0, {"x": 1000, "y": 1001, "w": 1002, "z": 1003}.get(val, 1100))
        if "unlocalized" in region:
            return (1, region)
        match = re.search(r'[_-](\d+)$', region)
        if match:
            return (2, int(match.group(1)))
        return (3, region)

    @staticmethod
    def weighted_stats(d):
        total = d["length"].sum()
        mean = (d["coverage"] * d["length"]).sum() / total
        d_sorted = d.sort_values("coverage")
        d_sorted["cumlen"] = d_sorted["length"].cumsum()
        median_row = d_sorted[d_sorted["cumlen"] >= total / 2].iloc[0]
        return {
            "mean": mean,
            "median": median_row["coverage"],
            "min": d["coverage"].min(),
            "max": d["coverage"].max(),
            "sequenced_bp": int((d["coverage"] * d["length"]).sum()),
            "total_bp": int(d["length"].sum()),
            "breadth_ge_1": int(d.loc[d["coverage"] >= 1, "length"].sum()),
            "breadth_ge_2": int(d.loc[d["coverage"] >= 2, "length"].sum()),
            "breadth_ge_5": int(d.loc[d["coverage"] >= 5, "length"].sum()),
            "breadth_ge_10": int(d.loc[d["coverage"] >= 10, "length"].sum())
        }

    @staticmethod
    def coverage_groups(d):
        total = d["length"].sum()

        def sum_bp(mask): return d.loc[mask, "length"].sum()

        return {
            "0× Coverage": (sum_bp(d["coverage"] == 0), total),
            "1–5× Coverage": (sum_bp((d["coverage"] >= 1) & (d["coverage"] <= 5)), total),
            "6–10× Coverage": (sum_bp((d["coverage"] >= 6) & (d["coverage"] <= 10)), total),
            ">10× Coverage": (sum_bp(d["coverage"] > 10), total)
        }

    @staticmethod
    def format_bp(bp):
        if bp >= 1_000_000:
            return f"{bp / 1_000_000:.2f} Mb"
        elif bp >= 1_000:
            return f"{bp / 1_000:.1f} kb"
        return f"{bp} bp"


    @staticmethod
    def make_cumulative_plot(df, max_x=100):
        cov = df["coverage"]
        weights = df["length"]
        bins = np.arange(0, min(max_x + 2, int(cov.max()) + 2))
        hist, edges = np.histogram(cov, bins=bins, weights=weights)
        total = weights.sum()
        cum = np.cumsum(hist[::-1])[::-1] / total * 100
        x_vals = edges[:-1]
        y_vals = cum

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=x_vals,
            y=y_vals,
            mode="lines",
            name="% of Genome",
            fill="tozeroy",
            line=dict(color="royalblue"),
            hovertemplate="Coverage ≥ %{x}×<br>% of Genome: %{y:.2f}%<extra></extra>"
        ))
        fig.update_layout(
            title="Cumulative Coverage Distribution",
            xaxis_title="Coverage ≥ x",
            yaxis_title="% of Genome",
            xaxis_range=[0, max_x],
            height=400,
            margin=dict(l=40, r=40, t=40, b=80),
            autosize=True
        )
        return fig

    @staticmethod
    def make_coverage_histogram(df, max_x=100):
        cov = df["coverage"]
        weights = df["length"]
        bins = np.arange(0, max_x + 2)
        hist, edges = np.histogram(cov, bins=bins, weights=weights)
        mids = bins[:-1]

        df_hist = pd.DataFrame({"Coverage": mids, "Bases": hist})
        fig = px.bar(df_hist, x="Coverage", y="Bases",
                     title="Coverage Distribution",
                     labels={"Coverage": "Coverage (×)", "Bases": "Bases"},
                     height=400)
        fig.update_layout( margin=dict(l=40, r=40, t=40, b=80), autosize=True)
        return fig

    @staticmethod
    def make_breadth_per_region_plot(df, min_coverage=1):
        regions = []
        breadths = []

        for region, sub in df.groupby("region"):
            total = sub["length"].sum()
            covered = sub.loc[sub["coverage"] >= min_coverage, "length"].sum()
            percent = 100 * covered / total if total else 0
            regions.append(region)
            breadths.append(percent)

        data = pd.DataFrame({"Region": regions, f"Breadth ≥{min_coverage}× (%)": breadths})
        sorted_regions = sorted(data["Region"], key=GenomeCoverage.region_sort_key)

        data["Region"] = pd.Categorical(data["Region"], categories=sorted_regions, ordered=True)
        data = data.sort_values("Region")

        fig = px.bar(
            data,
            x="Region",
            y=f"Breadth ≥{min_coverage}× (%)",
            title=f"Breadth ≥{min_coverage}× Coverage per Region",
            labels={"Region": "Chromosome / Contig"},
            height=400
        )
        fig.update_layout(xaxis_tickangle=45, yaxis=dict(range=[0, 100]), margin=dict(l=40, r=40, t=40, b=80), autosize=True)
        return fig

    @staticmethod
    def make_coverage_heatmap(df, bin_size=100000, zmax=None, colorscale="Inferno"):
        df = df.copy()
        df["bin"] = (df["start"] // bin_size) * bin_size

        grouped = df.groupby(["region", "bin"], as_index=False).apply(
            lambda g: pd.Series({
                "mean_coverage": (g["coverage"] * g["length"]).sum() / g["length"].sum()
            }), include_groups=False
        ).reset_index(drop=True)

        # Schritt 2: Länge pro Region → Min/Max Bin bestimmen
        region_ranges = df.groupby("region").agg({
            "start": "min",
            "end": "max"
        })
        region_ranges["bin_start"] = (region_ranges["start"] // bin_size) * bin_size
        region_ranges["bin_end"] = (region_ranges["end"] // bin_size) * bin_size

        # Schritt 3: Volles Raster (alle Bins pro Region) erzeugen
        full_bins = []
        for region, row in region_ranges.iterrows():
            for b in range(int(row["bin_start"]), int(row["bin_end"]) + bin_size, bin_size):
                full_bins.append((region, b))
        full_df = pd.DataFrame(full_bins, columns=["region", "bin"])

        # Schritt 4: Join mit berechneten Coverage-Werten
        merged = full_df.merge(grouped, on=["region", "bin"], how="left")
        merged["mean_coverage"] = merged["mean_coverage"].fillna(0)

        # Schritt 5: Pivot + NaN für nicht belegte Außenbereiche setzen
        pivot = merged.pivot(index="region", columns="bin", values="mean_coverage")

        # Maske: Bin liegt außerhalb des echten Bereichs → auf NaN setzen
        for region in pivot.index:
            bin_start = region_ranges.loc[region, "bin_start"]
            bin_end = region_ranges.loc[region, "bin_end"]
            for b in pivot.columns:
                if b < bin_start or b > bin_end:
                    pivot.at[region, b] = np.nan

        # Sortieren nach deiner Logik
        pivot = pivot.reindex(sorted(pivot.index, key=GenomeCoverage.region_sort_key))

        if zmax is None:
            zmax = np.nanquantile(pivot.values, 0.995)

        heatmap = go.Heatmap(
            z=pivot.values,
            x=pivot.columns,
            y=pivot.index,
            coloraxis="coloraxis",
            hoverongaps=False
        )

        fig = go.Figure(data=[heatmap])
        fig.update_layout(
            coloraxis=dict(
                cmin=0,
                cmax=zmax,
                colorscale=colorscale,
                colorbar=dict(title="Coverage (×)")
            ),
            title=f"Coverage Heatmap per chromosome/contig (bin size: {GenomeCoverage.format_bp(bin_size)})",
            xaxis_title="Genomic Position (binned)",
            yaxis_title="Region",
            height=700,
            margin=dict(l=40, r=40, t=60, b=80),
            plot_bgcolor='rgba(240,240,240,1)'
        )

        return fig

    @staticmethod
    def wrap_plot_in_accordion_with_slider(fig, div_id="myHeatmap", bin_size=100000,
                                           zmax_default=5, zmax_max=50) -> str:
        fig_json = json.dumps(fig, cls=PlotlyJSONEncoder)

        return f'''
      <div style="margin-top:10px;">
        <div style="margin-bottom: 10px;">
          <label for="{div_id}_slider"><b>Max. Coverage (zmax):</b> <span id="{div_id}_val">{zmax_default}</span>×</label><br>
          <input type="range" id="{div_id}_slider" min="1" max="{zmax_max}" value="{zmax_default}" step="1" style="width:300px;">
          <p class="plot-note" style="margin-top: 5px;">
            Adjust this slider to change the upper limit of the color scale. Lower values increase contrast in low-coverage regions, 
            while higher values reveal high-coverage peaks.
          </p>
        </div>
        <div id="{div_id}" style="width: 100%; height: 700px;"></div>
        <p class="plot-note">This heatmap provides a chromosome-by-position overview of coverage depth across the genome.
        Each cell represents a fixed-size genomic bin ({GenomeCoverage.format_bp(bin_size)}), colored by its average coverage.</p>
      </div>


    <script>
    document.addEventListener("DOMContentLoaded", function() {{
      const fig = {fig_json};

      const plotDiv = document.getElementById("{div_id}");
      const slider = document.getElementById("{div_id}_slider");
      const zval   = document.getElementById("{div_id}_val");

      Plotly.newPlot(plotDiv, fig.data, fig.layout, {{responsive: true, scrollZoom: true}});

      slider.addEventListener("input", function() {{
        const newZmax = parseFloat(slider.value);
        zval.textContent = newZmax;
        Plotly.relayout(plotDiv, {{"coloraxis.cmax": newZmax}});
      }});
    }});
    </script>
    '''

    @staticmethod
    def write_coverage_summary_table(f, sequenced_bp, total_bp):
        percent = 100 * sequenced_bp / total_bp if total_bp else 0
        f.write("<div class='stat-plot'>\n")
        f.write("<div class='stat-block-box'>\n")
        f.write("<table class='stat-table'><caption>Sequencing output</caption>\n")
        f.write("<thead><tr><th class='table-header'>Metric</th><th class='table-header'>Value</th></tr></thead>\n")
        f.write("<tbody>\n")
        f.write(f"<tr><td>Total Sequenced Bases</td><td>{GenomeCoverage.format_bp(sequenced_bp)}</td></tr>\n")
        f.write(f"<tr><td>Reference Genome Size</td><td>{GenomeCoverage.format_bp(total_bp)}</td></tr>\n")
        #f.write(f"<tr><td>Percentage Covered</td><td>{percent:.1f}%</td></tr>\n")
        f.write("</tbody></table>\n")
        f.write(
            "<p class='plot-note'>Sequencing yield compared to reference genome length.</p>\n")
        f.write("</div></div>\n")

    @staticmethod
    def write_summary_block(f, stats, cov_stats, show_seq_output=True):
        f.write("<div class='stat-plot-row'>\n")

        # Block: Coverage-Statistiken
        f.write("<div class='stat-plot'>\n")
        f.write("<div class='stat-block-box'>\n")
        f.write("<table class='stat-table'><caption>Coverage Statistics</caption>\n")
        f.write("<thead><tr><th class='table-header'>Metric</th><th class='table-header'>Value</th></tr></thead>\n")
        f.write("<tbody>\n")
        for key, val in stats.items():
            if key in ("sequenced_bp", "total_bp", "breadth_ge_1", "breadth_ge_2", "breadth_ge_5", "breadth_ge_10"):
                continue  # separat ausgegeben
            f.write(f"<tr><td>{key.capitalize()}</td><td>{val:.2f}×</td></tr>\n")
        f.write("</tbody></table>\n")
        f.write(
            "<p class='plot-note'>Statistical overview of how coverage is distributed across all mapped positions.</p>\n")
        f.write("</div></div>\n")

        # Block: Coverage-Gruppen (z.B. 0–5×)
        f.write("<div class='stat-plot'>\n")
        f.write("<div class='stat-block-box'>\n")
        f.write("<table class='stat-table'><caption>Coverage Categories (bp / %)</caption>\n")
        f.write("<thead><tr><th class='table-header'>Category</th><th class='table-header'>Value</th></tr></thead>\n")
        f.write("<tbody>\n")
        total_length = None
        for label, (count, total) in cov_stats.items():
            if total_length is None:
                total_length = total
            percent = 100 * count / total if total else 0
            f.write(f"<tr><td>{label}</td><td>{GenomeCoverage.format_bp(count)} ({percent:.1f}%)</td></tr>\n")
        f.write("</tbody></table>\n")
        f.write(
            "<p class='plot-note'>How much of the genome is covered at different ranges. Useful to assess how much is poorly covered or not covered at all.</p>\n")
        f.write("</div></div>\n")

        # Block: Breadth-Werte separat
        f.write("<div class='stat-plot'>\n")
        f.write("<div class='stat-block-box'>\n")
        f.write("<table class='stat-table'><caption>Breadth of Coverage (bp / %)</caption>\n")
        f.write(
            "<thead><tr><th class='table-header'>Threshold</th><th class='table-header'>Covered Bases</th></tr></thead>\n")
        f.write("<tbody>\n")
        breadth_keys = [("breadth_ge_1", "1×"), ("breadth_ge_2", "2×"), ("breadth_ge_5", "5×"),
                        ("breadth_ge_10", "10×")]
        for key, label in breadth_keys:
            if key in stats and total_length:
                covered = stats[key]
                percent = 100 * covered / total_length
                f.write(
                    f"<tr><td>Breadth ≥ {label}</td><td>{GenomeCoverage.format_bp(covered)} ({percent:.1f}%)</td></tr>\n")
        f.write("</tbody></table>\n")
        f.write(
            "<p class='plot-note'>Breadth indicates how much of the genome is covered at or above a given coverage depth. Higher breadth indicates more uniform and complete coverage.</p>\n")
        f.write("</div></div>\n")

        f.write("</div>\n")

    def render_html_report(self):
        print("Generating HTML report...")
        total_stats = self.weighted_stats(self.df)
        coverage_stats = self.coverage_groups(self.df)
        fig_cum = self.make_cumulative_plot(self.df, max_x=20)
        fig_hist = self.make_coverage_histogram(self.df, max_x=20)

        with open(self.output_html, "w") as f:
            f.write("<html><head><meta charset='utf-8'>\n")
            #f.write("<script src='https://cdn.plot.ly/plotly-2.27.0.min.js'></script>\n")
            self.write_style(f)
            f.write("</head><body>")
            f.write(f"<h1 class='caption'>Coverage Report – {self.sample_name} - {self.species}</h1>\n")
            f.write("<p style='margin-bottom: 20px;'>"
                    "  This report provides a detailed overview of genome coverage based on sequencing data from a field-collected sample."
                    "  It is intended to help assess the overall quality of the sample by evaluating how well and how evenly the genome is covered."
                    "  The included statistics and plots allow you to determine whether the sequencing data is sufficient to support downstream applications such as genome assembly or taxonomic profiling."
                    "  This overview is especially useful for making quick decisions in field settings — for example, to judge whether a sample warrants further processing or whether additional sequencing may be required."
                    "</p>")
            ### GENOME WIDE SUMMARY ####
            f.write("<h3>Genome-wide Summary</h3>\n")
            self.write_summary_block(f, total_stats, coverage_stats)

            ### MAPPED READS ####
            f.write("<h3>Mapped Reads</h3>")
            mapped_reads = self.read_counts_per_species.get(self.species, 0)
            f.write(f"<p>For <b>{self.species}</b>, a total of <b>{mapped_reads:,}</b> primary mapped reads were detected "
                    f"based on the reference alignment. This includes only reads that are mapped, non-secondary, and non-supplementary.</p>\n")
            # Sequencing output summary (verlagert von summary_block)
            GenomeCoverage.write_coverage_summary_table(f, sequenced_bp=total_stats["sequenced_bp"],
                                                        total_bp=total_stats["total_bp"])

            ### BREATH OF COVERAGE ####
            fig_breadth_1 = self.make_breadth_per_region_plot(self.df, min_coverage=1)
            fig_breadth_2 = self.make_breadth_per_region_plot(self.df, min_coverage=2)
            fig_breadth_5 = self.make_breadth_per_region_plot(self.df, min_coverage=5)
            fig_breadth_10 = self.make_breadth_per_region_plot(self.df, min_coverage=10)

            f.write("<h3>Breadth of Coverage per Region</h3>\n")
            f.write("<div>\n")
            f.write("<div class='tab-buttons'>\n")
            f.write("<button onclick='showBreadthTab(\"1\")'>≥1×</button>\n")
            f.write("<button onclick='showBreadthTab(\"2\")'>≥2×</button>\n")
            f.write("<button onclick='showBreadthTab(\"5\")'>≥5×</button>\n")
            f.write("<button onclick='showBreadthTab(\"10\")'>≥10×</button>\n")
            f.write("</div>\n")

            # Tab: ≥1×
            f.write("<div class='plot-tab' id='breadth-tab-1' style=' height:400px; display: block;'>\n")
            f.write(pio.to_html(fig_breadth_1, include_plotlyjs=True, full_html=False,
                                config={'responsive': True, 'scrollZoom': True}))
            f.write("</div>\n")

            # Tab: ≥2×
            f.write("<div class='plot-tab' id='breadth-tab-2' style=' height:400px; display: none;'>\n")
            f.write(pio.to_html(fig_breadth_2, include_plotlyjs=False, full_html=False,
                                config={'responsive': True, 'scrollZoom': True}))
            f.write("</div>\n")

            # Tab: ≥5×
            f.write("<div class='plot-tab' id='breadth-tab-5' style=' height:400px; display: none;'>\n")
            f.write(pio.to_html(fig_breadth_5, include_plotlyjs=False, full_html=False,
                                config={'responsive': True, 'scrollZoom': True}))
            f.write("</div>\n")

            # Tab: ≥10×
            f.write("<div class='plot-tab' id='breadth-tab-10' style=' height:400px; display: none;'>\n")
            f.write(pio.to_html(fig_breadth_10, include_plotlyjs=False, full_html=False,
                                config={'responsive': True, 'scrollZoom': True}))
            f.write("</div>\n")

            f.write(
                "<p class='plot-note'>This plot shows the percentage of each region that is covered by at least one read. "
                "Regions with high breadth are likely to be informative for downstream analyses such as (reference-guided) assembly or alignment-based profiling.<br>"
                "<i>The breadth of coverage refers to the percentage of genome bases sequenced at or above a given sequencing depth.</i></p>")

            # JavaScript zum Umschalten
            f.write("""
            <script>
            function showBreadthTab(thresh) {
              ['1', '2', '5', '10'].forEach(t => {
                const tab = document.getElementById('breadth-tab-' + t);
                if (tab) {
                  const visible = (t === thresh);
                  tab.style.display = visible ? 'block' : 'none';
                  if (visible) {
                    // Verzögertes Resize für Layout-Stabilität
                    setTimeout(() => {
                      tab.querySelectorAll('.js-plotly-plot').forEach(plot => Plotly.Plots.resize(plot));
                    }, 100);
                  }
                }
              });
            }
            </script>
            """)


            ### COVERAGE DISTRIBUTION ####
            f.write("<h3>Coverage distribution</h3>\n")
            f.write("<div>\n")
            f.write("<div class='tab-buttons'>\n")
            f.write("<button  onclick='showTab(\"cum\")'>Cumulative</button>\n")
            f.write("<button  onclick='showTab(\"hist\")'>Histogram</button>\n")
            f.write("</div>\n")
            f.write("<div class='plot-tab' id='plot-hist' style='display:none;'>\n")
            ### add plotlyjs in first plot
            f.write(pio.to_html(fig_hist, include_plotlyjs=False, full_html=False, config={'responsive': True, 'scrollZoom': True}))
            f.write(
                "<p class='plot-note'>This histogram shows how many bases in the genome are covered at each coverage level (e.g., 0×, 1×, 2×, ...). It helps identify the overall shape of the coverage profile: peaks around the expected depth indicate good sequencing quality, while a high fraction of bases with low or zero coverage may suggest dropouts, technical issues, or poorly mappable regions.</p>\n")
            f.write("</div>\n<div class='plot-tab' id='plot-cum' style='display:block;'>\n")
            f.write(pio.to_html(fig_cum, include_plotlyjs=False, full_html=False, config={'responsive': True, 'scrollZoom': True}))
            f.write(
                "<p class='plot-note'>This cumulative plot shows the percentage of bases in the genome that are covered at or above a given coverage threshold (e.g., 10×). It provides a quick overview of how much of the genome meets typical analysis requirements. For example, if 95% of bases have ≥10× coverage, variant calling is likely reliable for most regions.</p>\n")
            f.write("</div>\n</div>\n</div>\n")


            ### HEATMAP ####
            f.write("<h3>Coverage heatmap</h3>\n")
            heatmap_bin_size = 100000
            heatmap_fig = GenomeCoverage.make_coverage_heatmap(self.df, bin_size=heatmap_bin_size, zmax=5)
            html = GenomeCoverage.wrap_plot_in_accordion_with_slider(
                heatmap_fig, bin_size=heatmap_bin_size
            )
            f.write(html)

            regions_to_plot = self.regions_sorted
            # Begrenzung für Lyrurus tetrix
            limit_lyrurus_tetrix = 30
            if self.species == "Lyrurus tetrix" and len(regions_to_plot) > limit_lyrurus_tetrix:
                print(
                    f"⚠️  Warning: Showing only the {limit_lyrurus_tetrix} longest contigs for Lyrurus tetrix due to large genome size ({len(self.regions_sorted)} total regions).")

                # Länge pro Region berechnen und Top N auswählen
                region_lengths = (
                    self.df.groupby("region")["length"].sum()
                    .sort_values(ascending=False)
                    .head(30)
                    .index.tolist()
                )
                regions_to_plot = sorted(region_lengths, key=self.region_sort_key)
                f.write(f"<div class='plot-note' style='color: red; font-weight: bold;'>"
                        f"Warning: Only the {limit_lyrurus_tetrix} longest contigs are shown for <i>Lyrurus tetrix</i> "
                        f"due to large genome size (total: {len(self.regions_sorted)} contigs).</div>\n")
            # Navigation
            f.write("<h2 class='chrom-wise-header'>Chromosome-wise Coverage Analysis</h2>\n")
            f.write("<nav>\n<button onclick='prevPage()'>⟵ Previous</button>\n")
            f.write("<select id='pageSelect' onchange='goToPage(this.value)'>\n")
            for i, region in enumerate(regions_to_plot):
                f.write(f"<option value='{i}'>{region}</option>\n")
            f.write("</select>\n<button onclick='nextPage()'>Next ⟶</button>\n</nav>\n")
            f.write("<hr>")
            for i, region in enumerate(regions_to_plot):
                print(f"Processing {region}  ({i+1}/{len(regions_to_plot)})")
                sub = self.df[self.df["region"] == region].copy()
                sub["bin"] = (sub["start"] // self.bin_size) * self.bin_size
                binned = sub.groupby("bin", as_index=False).apply(
                    lambda g: pd.Series({
                        "coverage": (g["coverage"] * g["length"]).sum() / g["length"].sum()
                    }), include_groups=False
                ).reset_index()

                fig = px.bar(
                    binned,
                    x="bin",
                    y="coverage",
                    labels={"bin": "Position", "coverage": "Mean coverage"},
                    height=300,
                    opacity=1.0
                )
                fig.update_traces(marker_line_width=0)
                fig.update_layout(margin=dict(l=40, r=40, t=40, b=30))

                stats = self.weighted_stats(sub)
                cov_stats = self.coverage_groups(sub)

                f.write(f"<div class='page' id='page{i}'>\n")
                f.write(f"<h2>{region}</h2>\n")
                self.write_summary_block(f, stats, cov_stats, show_seq_output=False)
                f.write(f"<h3>Coverage plot ({self.bin_size // 1000} kb bins)</h3>")
                f.write(
                    f"<div class='plot'>{pio.to_html(fig, include_plotlyjs=False, full_html=False, config={'responsive': True, 'scrollZoom': True})}</div>\n")
                f.write(
                    f"<p class='plot-note'>This plot shows the average coverage across <b>{region}</b>, "
                    f"calculated in fixed-size bins ({GenomeCoverage.format_bp(self.bin_size)}). "
                    f"It provides a quick visual overview of how evenly the genome is covered along this region. "
                    f"This helps assess whether the sample delivers sufficient and consistent coverage — an important criterion for evaluating sample quality, "
                     f"especially in field settings or when working with limited or degraded input material. "
                    f"Note that due to binning, local fluctuations at single-base resolution may be smoothed or hidden. "
                    f"The plot is intended to highlight broad patterns, not per-base detail.</p>\n")
                f.write("</div>\n")
            f.write(f"""
                    <script>
                    let currentPage = 0;
                    const totalPages = {len(self.regions_sorted)};
                    const pages = document.querySelectorAll(".page");
                
                    function showPage(i) {{
                        pages.forEach(p => p.classList.remove("active"));
                        if (pages[i]) {{
                            pages[i].classList.add("active");
                            document.getElementById("pageSelect").value = i;
                            currentPage = i;
                            pages[i].querySelectorAll(".js-plotly-plot").forEach(p => Plotly.Plots.resize(p));
                        }}
                    }}
                    function showTab(tabName) {{
                      ['hist', 'cum'].forEach(name => {{
                        const el = document.getElementById(`plot-${{name}}`);
                        if (el) {{
                          const visible = (name === tabName);
                          el.style.display = visible ? 'block' : 'none';
                          if (visible) {{
                            setTimeout(() => {{
                              el.querySelectorAll('.js-plotly-plot').forEach(p => Plotly.Plots.resize(p));
                            }}, 100);
                          }}
                        }}
                      }});
                    }}
                    function prevPage() {{ if (currentPage > 0) showPage(currentPage - 1); }}
                    function nextPage() {{ if (currentPage < totalPages - 1) showPage(currentPage + 1); }}
                    function goToPage(i) {{ showPage(parseInt(i)); }}
                    document.addEventListener("DOMContentLoaded", () => showPage(0));
                    </script>
                    """)
            f.write("""
            <footer style="
              background: #333333;
              color: #eee;
              text-align: center;
              padding: 12px 0;
              margin-top: 60px;
              font-size: 13px;
              font-family: sans-serif;
              width: 100%;
            ">
              Created 2025 by Dr. Alexander Graf – feedback & issues: 
                <a href="mailto:alexander.graf@lmu.de" style="color: #80bfff; text-decoration: none;">alexander.graf@lmu.de</a>
            </footer>
            </body></html>
            """)
            print(f"HTML report generated in: {self.output_html}")

    @staticmethod
    def write_style(f):
        f.write("""
           <style>
           body { font-family: sans-serif; margin: 30px; }
           .plot { width: 100%; min-height: 320px; }
           .page { display: none; }
           .page.active { display: block; }
           nav { margin: 20px 0; }
           h2 { margin-top: 40px; }
           h3 { background: #c5c5c5; padding: 10px; }
           .caption { background: #6dadff; padding: 20px; }
           .chrom-wise-header{
                background: #2f2f2f;
                color: white;
                padding: 20px;
           }

           .stat-plot-row {
               display: flex;
               flex-wrap: nowrap;
               gap: 40px;
               margin: 20px 0;
               justify-content: flex-start;
               align-items: flex-start;
           }

           .stat-plot {
                flex: 1 1 1200px;
                max-width: 1200px;
                min-width: 300px;
            }

           .stat-table {
               border-collapse: collapse;
               font-size: 14px;
               min-width: 340px;
           }
           .table-header{
                background:#6dadff;
           }
           .stat-table caption {
               font-weight: bold;
               text-align: left;
               margin-bottom: 6px;
               font-size: 15px;
           }
           .stat-table th, .stat-table td {
               padding: 6px 10px;
               border: 1px solid #ccc;
           }
           .stat-table td:first-child { text-align: left; }
           .stat-table td:last-child { text-align: right; font-family: monospace; }

           .tab-buttons { margin-bottom: 10px; }
           
           .plot-tab { width: 100%;  min-height: 320px; }
           .plot-tab .js-plotly-plot { width: 100% !important; height: 400px !important; }
           .plot-note {
                font-size: 13px;
                color: #555;
                margin-top: 0px;
                line-height: 1.5;
            }
            button {
              background-color: #1e88e5;  /* kräftiges Blau */
              color: white;
              border: none;
              border-radius: 5px;
              padding: 6px 14px;
              font-size: 14px;
              cursor: pointer;
              margin: 5px;
              font-weight: bold;
              transition: background-color 0.2s ease;
            }
            
            button:hover {
              background-color: #1565c0;  /* dunkleres Blau beim Hover */
            }
            select {
              background-color: #f0f0f0;
              border: 1px solid #ccc;
              border-radius: 5px;
              padding: 6px 10px;
              font-size: 14px;
              color: #333;
              font-family: sans-serif;
              cursor: pointer;
              margin: 5px;
              transition: border-color 0.2s ease, box-shadow 0.2s ease;
            }
            
            select:focus {
              border-color: #1e88e5;
              box-shadow: 0 0 4px rgba(30, 136, 229, 0.4);
              outline: none;
            }

           </style>
           """)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Coverage report with disjoint coverage stats and barplots.")
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--output_dir", default="coverage_report_html")
    parser.add_argument("--species", required=True,  choices=["Lagopus muta", "Lepus timidus", "Lyrurus tetrix"])
    parser.add_argument("--bin_size", type=int, default=10000)
    args = parser.parse_args()

    report = GenomeCoverage(
        mapping_file=args.mapping,
        sample_name=args.sample,
        bed_file=args.bed,
        bam_file= args.bam,
        species=args.species,
        bin_size=args.bin_size,
        output_dir=args.output_dir
    )
    report.render_html_report()