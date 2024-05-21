#!/usr/bin/env python3

import os
import sys
import pandas
import numpy
import plotly.express as px
import plotly.graph_objects as go
import itertools
import functools
import scipy

# ---- Load inputs ----

input_table_filt = sys.argv[1]
input_table_full = sys.argv[2]
design_tab = sys.argv[3]
group_A = sys.argv[4]
group_B = sys.argv[5]
output_dir = sys.argv[6]
min_comparison_pct = int(sys.argv[7])
min_dpsi = int(sys.argv[8])

# ---- Utility functions ----

FIG_TEMPLATE = dict(
    layout=go.Layout(
        autosize=False,
        title_font=dict(family="Arial", size=16),
        font=dict(family="Arial", size=16, color="#1e2125"),
        xaxis=dict(
            showgrid=False,
            ticks="outside",
            ticklen=10,
            tickwidth=2,
            showline=True,
            linecolor="black",
            tickcolor="black",
            mirror=True,
            linewidth=2,
            zeroline=False,
        ),
        yaxis=dict(
            showgrid=False,
            ticks="outside",
            ticklen=5,
            tickwidth=2,
            showline=True,
            linecolor="black",
            tickcolor="black",
            mirror=True,
            linewidth=2,
            zeroline=False,
            nticks=10,
        ),
        legend=dict(
            font=dict(size=14),
            x=0.5,
            y=1.15,
            orientation="h",
            xanchor="auto",
            yanchor="auto",
        ),
    )
)


title_dict = dict(
    EX="Exon Inclusion Levels",
    INT="Intron Retention Levels",
    ALT="Alternative Splice Site Usage",
)

axes_dict = dict(
    EX="PSI<br>[Percent Spliced In, %]",
    INT="PRI<br>[Percent Retained In, %]",
    ALT="PSU<br>[Percent SpliceSite Usage, %]",
)


def plot_boxplots(
    comparisons_data: pandas.DataFrame,
    filtered_data: pandas.DataFrame,
    event_type: str,
    group_a: str,
    group_b: str,
    min_pct_comparisons: int,
    min_dpsi: float,
):
    """

    Plot scatter of splicing quantification levels for a given event type.

    Parameters
    ----------
    comparisons_data : pandas.DataFrame
        Dataframe of pairwise comparisons.
    filtered_data : pandas.DataFrame
        Filtered dataframe pairwise comparisons filtered, in long format (melt).
    event_type : str
        Event type for data selection.
    group_a : str
        Control condition.
    group_b : str
        Contrast condition.
    min_pct_comparisons: int
        Minimum percentage of valid comparisons per event across all comparisons.

    Returns
    -------
    plotly.graph_objs.Figure
        Scatter plot of splicing quantification levels.

    """

    medians = px.line(
        comparisons_data.sort_values("MEDIAN_DPSI")[["EVENT", "MEDIAN_DPSI"]],
        x="EVENT",
        y="MEDIAN_DPSI",
        title="",
        color_discrete_sequence=["darkred"],
    )
    n_events = filtered_data["EVENT"].nunique()
    fig_boxplots = px.box(
        filtered_data,
        x="EVENT",
        y="dPSI",
        color="GROUP",
        notched=True,
        width=1000,
        height=600,
        color_discrete_map={
            "Down": "cornflowerblue",
            "NonDiff": "lightgrey",
            "Up": "darkorange",
        },
        labels={
            "EVENT": f"Regulated {event_type}",
            "dPSI": f"\u0394PSI<br>{group_b} vs {group_a}",
            "GROUP": "Regulation:",
        },
        title=f"{n_events} regulated events ({event_type}) in at least {min_pct_comparisons}% of comparisons <br>and min. median \u0394PSI +/- {min_dpsi}",
    )
    fig_boxplots.add_trace(medians.data[0])
    fig_boxplots.update_layout(
        template=FIG_TEMPLATE,
        xaxis=dict(showticklabels=False, ticks=""),
        margin=dict(pad=0, t=150, l=100),
    )
    fig_boxplots.update_traces(marker=dict(opacity=0.75, size=3))

    return fig_boxplots


def plot_scatter(
    medians_data: pandas.DataFrame,
    diff_events: list,
    event_type: str,
    group_a: str,
    group_b: str,
    color_unreg: str,
    color_reg: str,
):
    """

    Parameters
    ----------
    medians_data: pandas.DataFrame
            Dataframe with EventID and per-group median PSI values.
    diff_events : list
            List of regulated events.
    event_type : str
            Type of events: EX - exons, INT - introns, ALT - alternative splice sites.
    group_a : str
            Control condition.
    group_b : str
            Contrast condition.
    color_unreg : str
            Color for unregulated events.
    color_reg : str
            Color for regulated events.

    Returns
    -------
        plotly.graph_objs.Figure
            Scatter plot of splicing quantification levels.
    """
    events_medians = (
        medians_data.set_index("EVENT").filter(regex=event_type, axis=0).reset_index()
    )
    events_medians["DIFF"] = events_medians["EVENT"].apply(
        lambda x: "Yes" if x in diff_events else "No"
    )

    fig_scatter = px.scatter(
        events_medians,
        x=group_a,
        y=group_b,
        color="DIFF",
        color_discrete_sequence=[color_unreg, color_reg],
        title=event_type,
        labels={"DIFF": "Regulated"},
    )
    fig_scatter.update_layout(
        title={
            "text": title_dict[event_type],
            "y": 0.91,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
        template=FIG_TEMPLATE,
        width=700,
        height=700,
        xaxis=dict(title=f"{group_a} {axes_dict[event_type]}"),
        yaxis=dict(title=f"{group_b} {axes_dict[event_type]}"),
        showlegend=True,
        margin=dict(pad=0, l=100, r=50, b=100, t=10),
    )
    fig_scatter.update_traces(
        marker=dict(size=9, line=dict(width=1, color=True)),
        selector=dict(mode="markers"),
        opacity=0.6,
    )
    fig_scatter.add_shape(
        type="line", x0=0.5, y0=0, x1=100, y1=100, line=dict(color="black", width=2)
    )

    return fig_scatter


# ---- SETUP ----
output_dir = output_dir + "/" + group_B + "_vs_" + group_A
os.makedirs(output_dir, exist_ok=True)
os.chdir(output_dir)

os.makedirs("figures", exist_ok=True)
os.makedirs("tables", exist_ok=True)

# Remove cryptic and constitutive events
stats = (
    pandas.read_csv(input_table_filt, sep="\t")
    .query("TYPE != 'HIGH_PSI' & TYPE != 'LOW_PSI'")
    .query("N_AS >= 5")
)
events = stats[["EVENT", "GENE", "COORD", "LENGTH"]]

psi_table = pandas.read_csv(input_table_full, sep="\t").drop(
    ["GENE", "LENGTH", "COORD", "FullCO", "COMPLEX"], axis=1
)

psi_table = psi_table.drop(psi_table.filter(regex="-Q").columns.values, axis=1)
psi_table = psi_table[psi_table["EVENT"].isin(stats["EVENT"].tolist())]
psi_table.to_csv("tables/PSI_TABLE.tab", sep="\t", index=False)

# Select samples for pair-wise comparisons
groups = pandas.read_csv(design_tab, sep="\t")

SAMPLES_A = groups[groups["GroupName"].str.contains(f"{group_A}")]["SampleName"].tolist()
SAMPLES_B = groups[groups["GroupName"].str.contains(f"{group_B}")]["SampleName"].tolist()
SAMPLES = SAMPLES_A + SAMPLES_B

# Generate a list of pair-wise comparisons
COMPARISONS = list(itertools.product(SAMPLES_A, SAMPLES_B))

# Set up parameters for the pair-wise comparisons
N_COMPARISONS = len(COMPARISONS)
Q_FILTER = "Q1 < 0 & Q3 < 0 | Q1 > 0 & Q3 > 0"

# Modify colors of regulated events for scatter plots
REGULATED_COLORS = dict(EX="teal", INT="blue", ALT="orange")

# Select samples from the vast-tools table
selected_data = psi_table.set_index("EVENT")[SAMPLES].reset_index()

# Calculate median PSI per group
medians = selected_data.set_index("EVENT")
medians[f"{group_A}"] = medians[SAMPLES_A].median(axis=1).round(2)
medians[f"{group_B}"] = medians[SAMPLES_B].median(axis=1).round(2)
medians_table = medians.reset_index()[["EVENT", f"{group_A}", f"{group_B}"]]

# ---- COMPARISONS ----
# Perform pairwise comparisons
PAIRWISE_COMPARISONS = []
for comparison in COMPARISONS:
    control = comparison[0]
    contrast = comparison[1]
    df_temp = selected_data.copy()[["EVENT", control, contrast]]
    df_temp[f"{contrast}-vs-{control}"] = round(df_temp[contrast] - df_temp[control], 2)
    df_temp.drop([control, contrast], axis=1, inplace=True)
    PAIRWISE_COMPARISONS.append(df_temp)

# Combine pairwise comparison data
comparisons_data = functools.reduce(
    lambda left, right: pandas.merge(left, right, on="EVENT"), PAIRWISE_COMPARISONS
)

# Get stats for the pair-wise comparisons
comparisons_data["N_VALID"] = (
    comparisons_data.filter(regex="vs")
    .apply(lambda x: abs(x) >= min_dpsi, axis=1)
    .sum(axis=1)
)
comparisons_data["PCT_VALID"] = round(
    (comparisons_data["N_VALID"] * 100) / N_COMPARISONS, 2
)
comparisons_data["MEDIAN_DPSI"] = (
    comparisons_data.filter(regex="vs")
    .apply(lambda x: numpy.quantile(x, q=0.5), axis=1)
    .round(2)
)
comparisons_data["Q1"] = (
    comparisons_data.filter(regex="vs")
    .apply(lambda x: numpy.quantile(x, q=0.25), axis=1)
    .round(2)
)
comparisons_data["Q3"] = (
    comparisons_data.filter(regex="vs")
    .apply(lambda x: numpy.quantile(x, q=0.75), axis=1)
    .round(2)
)
comparisons_data["IQR"] = (
    comparisons_data.filter(regex="vs")
    .apply(lambda x: scipy.stats.iqr(x), axis=1)
    .round(2)
)

comparisons_data_final = comparisons_data.merge(events, on="EVENT", how="left").merge(
    medians_table, on="EVENT", how="left"
)

comparisons_data_final.to_csv(
    f"tables/PAIRWISE_COMPARISONS-{group_B}-{group_A}_dPSI-{min_dpsi}-VLOW.tab",
    sep="\t",
    index=False,
)

# Iterate over event types and generate data and figures
for EVENT_TYPE in ["EX", "INT", "ALT"]:

    # Filter for N% of valid comparisons
    filtered_comparisons_data = (
        comparisons_data_final.query(f"PCT_VALID >= {min_comparison_pct}")
        .query(Q_FILTER)
        .set_index("EVENT")
        .filter(regex=EVENT_TYPE, axis=0)
        .reset_index()
    )

    filtered_comparisons_data["GROUP"] = filtered_comparisons_data["MEDIAN_DPSI"].apply(
        lambda x: "Up" if x >= min_dpsi else "Down" if x <= -min_dpsi else "NonDiff"
    )

    # Get info for filtered comparison data - save
    filtered_comparisons_data_stats = filtered_comparisons_data.drop(
        filtered_comparisons_data.filter(regex="vs", axis=1).columns, axis=1
    ).sort_values("MEDIAN_DPSI")

    # Save stats for filtered comparisons
    filtered_comparisons_data_stats.to_csv(
        f"tables/Diff{EVENT_TYPE}_MinPct-{min_comparison_pct}_MindPSI-{min_dpsi}.tab",
        sep="\t",
        index=False,
    )

    # Sort events by median dPSI
    sorter = filtered_comparisons_data.sort_values("MEDIAN_DPSI")["EVENT"].tolist()

    # Melt filter data for visualisation
    filtered_data_melt = (
        filtered_comparisons_data.set_index(["EVENT", "GROUP"])
        .filter(regex="vs")
        .reset_index()
        .melt(["EVENT", "GROUP"], value_name="dPSI", var_name="COMPARISON")
    )
    filtered_data_melt["EVENT"] = pandas.Categorical(
        filtered_data_melt["EVENT"], categories=sorter, ordered=True
    )
    filtered_data_melt.sort_values("EVENT", inplace=True)

    # List of regulated events
    events_diff = set(filtered_data_melt["EVENT"].tolist())

    dpsi_boxplots = plot_boxplots(
        filtered_comparisons_data,
        filtered_data_melt,
        event_type=EVENT_TYPE,
        group_a=group_A,
        group_b=group_B,
        min_pct_comparisons=min_comparison_pct,
        min_dpsi=min_dpsi,
    )

    psi_scatter = plot_scatter(
        medians_data=medians_table,
        diff_events=events_diff,
        event_type=EVENT_TYPE,
        group_a=group_A,
        group_b=group_B,
        color_reg=REGULATED_COLORS[EVENT_TYPE],
        color_unreg="lightgray",
    )

    dpsi_boxplots.write_image(
        file=f"figures/Diff{EVENT_TYPE}_PW-MinPct-{min_comparison_pct}_MindPSI-{min_dpsi}.pdf",
        format="pdf",
    )

    psi_scatter.write_image(
        file=f"figures/ScatterPlot_{EVENT_TYPE}_PW-MinPct-{min_comparison_pct}_MindPSI-{min_dpsi}.pdf",
        format="pdf",
    )

# Save a full table with all filtered comparisons across all event types
full_filtered_df = []
for EVENT_TYPE in ["EX", "INT", "ALT"]:
    df = pandas.read_csv(
        f"tables/Diff{EVENT_TYPE}_MinPct-{min_comparison_pct}_MindPSI-{min_dpsi}.tab",
        sep="\t",
    )
    df["EVENT_TYPE"] = EVENT_TYPE
    full_filtered_df.append(df)

full_filtered_df = pandas.concat(full_filtered_df, axis=0)
full_filtered_df.to_csv(
    f"tables/ALL_EVENTS_MinPct-{min_comparison_pct}_MindPSI-{min_dpsi}.tab",
    sep="\t",
    index=False,
)
