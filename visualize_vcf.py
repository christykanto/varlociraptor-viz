import altair as alt
import pandas as pd
import pysam
import re


def phred_to_prob(phred_value):
    """Convert PHRED score to probability"""
    if isinstance(phred_value, (tuple, list)):
        phred_value = phred_value[0]
    return 10 ** (-phred_value / 10)


def visualize_event_probabilities(record):
    """Visualize event probabilities from INFO column (PROB_* fields)"""
    prob_data = []

    for key, value in record.info.items():
        if key.startswith("PROB_"):
            event_name = key.replace("PROB_", "")

            probability = 0.0 if value == float("inf") else phred_to_prob(value)

            prob_data.append({"Event": event_name, "Probability": probability})

    df = pd.DataFrame(prob_data)

    chart = (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X("Event:N", title="Event Type"),
            y=alt.Y("Probability:Q", title="Probability"),
            tooltip=["Event", alt.Tooltip("Probability:Q", format=".6f")],
        )
        .properties(title="Event Probabilities", width=400, height=300)
    )

    return chart


def visualize_allele_frequency_distribution(record, sample_name):
    """Visualize allele frequency distribution (AFD field)"""
    sample = record.samples[sample_name]

    afd = sample["AFD"]
    af_ml = sample["AF"]

    if isinstance(af_ml, (tuple, list)):
        af_ml = af_ml[0]

    afd_entries = afd if isinstance(afd, (tuple, list)) else [afd]

    afd_data = []

    for entry in afd_entries:
        if isinstance(entry, str):
            for part in entry.split(","):
                if "=" in part:
                    freq, phred = part.split("=")
                    freq = float(freq)
                    prob = phred_to_prob(float(phred))

                    afd_data.append(
                        {
                            "Allele Frequency": freq,
                            "Probability": prob,
                            "Type": (
                                "ML Estimate"
                                if abs(freq - af_ml) < 0.001
                                else "Distribution"
                            ),
                        }
                    )

    df = pd.DataFrame(afd_data)

    base = (
        alt.Chart(df[df["Type"] == "Distribution"])
        .mark_circle(size=60, opacity=0.7)
        .encode(
            x="Allele Frequency:Q",
            y=alt.Y("Probability:Q", axis=None),
            color=alt.value("blue"),
            tooltip=["Allele Frequency", "Probability", "Type"],
        )
    )

    ml = (
        alt.Chart(df[df["Type"] == "ML Estimate"])
        .mark_circle(size=100)
        .encode(
            x="Allele Frequency:Q",
            y=alt.Y("Probability:Q", axis=None),
            color=alt.value("red"),
            tooltip=["Allele Frequency", "Probability", "Type"],
        )
    )

    return (
        base + ml
    ).properties(
        title="Allele Frequency Distribution (ML Estimate in Red)",
        width=500,
        height=300,
    ).configure_view(strokeWidth=0).configure_axis(grid=False)


def visualize_observations(record, sample_name):
    """Visualize observations from OBS field"""
    sample = record.samples[sample_name]
    obs = sample["OBS"]

    obs_string = obs[0] if isinstance(obs, (tuple, list)) else obs

    ref_observations = []
    alt_observations = []

    pattern = r"(\d+)([a-zA-Z]{2})(.{8})"
    matches = re.findall(pattern, obs_string)

    for idx, match in enumerate(matches):
        count = int(match[0])
        odds_code = match[1]
        rest = match[2]

        allele_type = odds_code[0].upper()

        obs_entry = {
            "obs_index": idx,
            "count": count,
            "Posterior Odds": odds_code[1],
            "Strand": rest[3],
            "Read Position": rest[5],
            "Orientation": rest[4],
            "Softclip": rest[6],
            "Indel": rest[7],
            "Edit Distance": int(rest[0]) if rest[0].isdigit() else 0,
        }

        if allele_type == "A":
            alt_observations.append(obs_entry)
        else:
            ref_observations.append(obs_entry)

    metrics = [
        "Posterior Odds",
        "Strand",
        "Read Position",
        "Orientation",
        "Softclip",
        "Indel",
        "Edit Distance",
    ]

    max_count = max(
        sum(o["count"] for o in ref_observations) or 0,
        sum(o["count"] for o in alt_observations) or 0,
    )

    def create_panel(observations, allele, show_legend=True):
        rows = []
        for obs in observations:
            for metric in metrics:
                rows.append(
                    {
                        "Metric": metric,
                        "Category": str(obs[metric]),
                        "Count": obs["count"],
                        "obs_index": obs["obs_index"],
                    }
                )

        df = pd.DataFrame(rows)

        # ---- Edit distance domain logic ----
        edit_df = df[df["Metric"] == "Edit Distance"]

        edit_domain = None
        if not edit_df.empty:
            values = edit_df["Category"].astype(int).unique()
            if len(values) == 1:
                k = int(values[0])
                edit_domain = [0, k] if k > 0 else [0, 1]

        base = alt.Chart(df).encode(
            x=alt.X("Metric:N", sort=metrics),
            y=alt.Y("Count:Q", stack="zero", scale=alt.Scale(domain=[0, max_count])),
            order="obs_index:Q",
            tooltip=["Metric", "Category", "Count"],
        )

        edit_layer = (
            base.transform_filter(alt.datum.Metric == "Edit Distance")
            .mark_bar(size=18)
            .encode(
                color=alt.Color(
                    "Category:Q",
                    scale=alt.Scale(
                domain=[0, edit_domain[1]] if edit_domain else None,
                range=["#ff0000", "#FFCCCC"],  
            ),
                    legend=(
                        alt.Legend(title="Edit distance", type="gradient")
                        if show_legend
                        else None
                    ),
                )
            )
        )

        other_layer = (
            base.transform_filter(alt.datum.Metric != "Edit Distance")
            .mark_bar(size=18)
            .encode(
                color=alt.Color(
                    "Category:N",
                    legend=alt.Legend(title="Category") if show_legend else None,
                )
            )
        )

        return (edit_layer + other_layer).properties(
            width=220,
            height=400,
            title=f"{allele} Allele Observations",
        )

    ref_chart = create_panel(ref_observations, "REF", show_legend=False)
    alt_chart = create_panel(alt_observations, "ALT", show_legend=True)

    return alt.hconcat(ref_chart, alt_chart, spacing=10).resolve_scale(
        y="shared"
    ).configure_view(strokeWidth=0).configure_axis(grid=False)


if __name__ == "__main__":
    vcf = pysam.VariantFile("examples/example.vcf")

    record = next(vcf)
    sample_name = list(record.samples.keys())[0]

    visualize_event_probabilities(record).save("event_probabilities.html")
    visualize_allele_frequency_distribution(record, sample_name).save(
        "allele_frequency_distribution.html"
    )
    visualize_observations(record, sample_name).save("observations.html")

    print("✓ All charts generated successfully!")

    vcf.close()
