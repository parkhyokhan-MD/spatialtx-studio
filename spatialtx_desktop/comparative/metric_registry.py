from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class MetricDefinition:
    internal_name: str
    display_name: str
    category: str
    unit: str
    scale_sensitive: bool = False
    normalization_denominator: str = ""
    plot_group: str = ""
    interpretation_priority: int = 3
    observational_only: bool = False
    deprecated: bool = False
    delta_name: str = ""

    def to_dict(self) -> dict:
        return asdict(self)

    @property
    def resolved_delta_name(self) -> str:
        return self.delta_name or f"delta_{self.internal_name}"


PROGRAM = "Program scores"
TRANSITION = "Gradient and transition"
GRAPH = "Graph and adjacency"
TOPOLOGY = "Topology and fragmentation"
SCALE = "Sample scale and QC"
CONTEXT = "Observational context"


# Ordering here is the single deterministic ordering used by exports and plots.
METRIC_REGISTRY: tuple[MetricDefinition, ...] = (
    MetricDefinition("C_mean", "C program mean", PROGRAM, "score", plot_group="program_compatibility", interpretation_priority=2, delta_name="delta_C"),
    MetricDefinition("S_mean", "S program mean", PROGRAM, "score", plot_group="program_compatibility", interpretation_priority=2, delta_name="delta_S"),
    MetricDefinition("R_mean", "R balance-field mean", PROGRAM, "score", plot_group="program_compatibility", interpretation_priority=2, delta_name="delta_R"),
    MetricDefinition("C_median", "C program median", PROGRAM, "score", plot_group="program", interpretation_priority=2, delta_name="delta_C_median"),
    MetricDefinition("S_median", "S program median", PROGRAM, "score", plot_group="program", interpretation_priority=2, delta_name="delta_S_median"),
    MetricDefinition("R_median", "R balance-field median", PROGRAM, "score", plot_group="program", interpretation_priority=2, delta_name="delta_R_median"),
    MetricDefinition("R_std", "R variability", PROGRAM, "score", interpretation_priority=2, delta_name="delta_R_variability"),
    MetricDefinition("gradient_mean", "Gradient mean", TRANSITION, "score", plot_group="transition", interpretation_priority=2),
    MetricDefinition("gradient_q90", "Gradient 90th percentile", TRANSITION, "score", plot_group="transition", interpretation_priority=2),
    MetricDefinition("localized_interface_fraction", "Localized interface-like fraction", TRANSITION, "fraction", plot_group="transition", interpretation_priority=1),
    MetricDefinition("diffuse_fraction", "Diffuse transition fraction", TRANSITION, "fraction", plot_group="transition", interpretation_priority=1),
    MetricDefinition("transition_burden_score", "Transition burden score", TRANSITION, "score", plot_group="transition", interpretation_priority=1),
    MetricDefinition("R_crossing_fraction", "R zero-crossing fraction", TRANSITION, "fraction", plot_group="transition", interpretation_priority=1),
    MetricDefinition("adj_same_fraction", "Same-side adjacency fraction", GRAPH, "fraction", plot_group="graph", interpretation_priority=2),
    MetricDefinition("adj_zero_fraction", "Near-zero adjacency fraction", GRAPH, "fraction", plot_group="graph", interpretation_priority=2),
    MetricDefinition("adj_opposite_fraction", "Opposite-side adjacency fraction", GRAPH, "fraction", plot_group="graph", interpretation_priority=2),
    MetricDefinition("n_diffuse_components", "Diffuse component count", TOPOLOGY, "count", True, "n_valid_spots", "topology_raw", 4),
    MetricDefinition("n_small_components", "Small diffuse component count", TOPOLOGY, "count", True, "n_valid_spots", "topology_raw", 4),
    MetricDefinition("n_interface_components", "Interface component count", TOPOLOGY, "count", True, "n_valid_spots", "topology_raw", 4),
    MetricDefinition("interface_fragmentation_index", "Interface fragmentation index", TOPOLOGY, "ratio", plot_group="topology_normalized", interpretation_priority=2),
    MetricDefinition("small_component_fraction", "Small-component fraction", TOPOLOGY, "fraction", plot_group="topology_normalized", interpretation_priority=1),
    MetricDefinition("largest_diffuse_component_ratio", "Largest diffuse-component ratio", TOPOLOGY, "ratio", plot_group="topology_normalized", interpretation_priority=1),
    MetricDefinition("diffuse_components_per_1000_valid_spots", "Diffuse components per 1,000 valid spots", TOPOLOGY, "per 1,000 spots", False, "n_valid_spots", "topology_normalized", 1),
    MetricDefinition("diffuse_components_per_1000_in_tissue_spots", "Diffuse components per 1,000 in-tissue spots", TOPOLOGY, "per 1,000 spots", False, "n_in_tissue_spots", "topology_normalized", 2),
    MetricDefinition("diffuse_components_per_tissue_component", "Diffuse components per tissue component", TOPOLOGY, "ratio", False, "tissue_component_count", "topology_normalized", 1),
    MetricDefinition("small_components_per_1000_valid_spots", "Small components per 1,000 valid spots", TOPOLOGY, "per 1,000 spots", False, "n_valid_spots", "topology_normalized", 1),
    MetricDefinition("transition_components_per_1000_transition_spots", "Diffuse components per 1,000 transition spots", TOPOLOGY, "per 1,000 spots", False, "n_transition_spots", "topology_normalized", 1),
    MetricDefinition("interface_segments_per_1000_valid_spots", "Interface components per 1,000 valid spots", TOPOLOGY, "per 1,000 spots", False, "n_valid_spots", "topology_normalized", 2),
    MetricDefinition("normalized_fragmentation_score", "Normalized transition-component density", TOPOLOGY, "per 1,000 spots", False, "n_valid_spots", "topology_normalized", 1),
    MetricDefinition("n_total_spots", "Total spots", SCALE, "spots", True, "", "sample_scale", 4),
    MetricDefinition("n_valid_spots", "Valid analysis spots", SCALE, "spots", True, "", "sample_scale", 3),
    MetricDefinition("n_in_tissue_spots", "In-tissue spots", SCALE, "spots", True, "", "sample_scale", 3),
    MetricDefinition("n_transition_spots", "Transition-candidate spots", SCALE, "spots", True, "", "sample_scale", 4),
    MetricDefinition("tissue_area_proxy", "Tissue area proxy (valid spots)", SCALE, "spot-count proxy", True, "", "sample_scale", 3),
    MetricDefinition("tissue_component_count", "Tissue graph components", SCALE, "components", True, "", "sample_scale", 3),
    MetricDefinition("spatial_extent_x", "Spatial extent X", SCALE, "coordinate units", True, "", "sample_scale", 3),
    MetricDefinition("spatial_extent_y", "Spatial extent Y", SCALE, "coordinate units", True, "", "sample_scale", 3),
    MetricDefinition("spatial_extent_area_proxy", "Spatial bounding-box area proxy", SCALE, "coordinate units squared", True, "", "sample_scale", 3),
    MetricDefinition("mean_spot_spacing", "Mean nearest-spot spacing", SCALE, "coordinate units", False, "", "sample_scale", 3),
    MetricDefinition("H_expr_mean", "H centered sample mean", CONTEXT, "centered score", False, "", "", 5, True, True, "delta_H_expr"),
    MetricDefinition("V_expr_mean", "V centered sample mean", CONTEXT, "centered score", False, "", "", 5, True, True, "delta_V_expr"),
    MetricDefinition("H_raw_mean", "H raw-scale mean", CONTEXT, "expression score", False, "", "hv", 3, True),
    MetricDefinition("H_raw_median", "H raw-scale median", CONTEXT, "expression score", False, "", "hv", 2, True),
    MetricDefinition("H_q75", "H 75th percentile", CONTEXT, "expression score", False, "", "hv", 3, True),
    MetricDefinition("H_q90", "H 90th percentile", CONTEXT, "expression score", False, "", "hv", 2, True),
    MetricDefinition("H_high_fraction", "H pooled-threshold high fraction", CONTEXT, "fraction", False, "pooled_reference_target_q90", "hv", 2, True),
    MetricDefinition("H_variance", "H variance", CONTEXT, "variance", False, "", "hv", 3, True),
    MetricDefinition("H_MAD", "H median absolute deviation", CONTEXT, "expression score", False, "", "hv", 3, True),
    MetricDefinition("H_transition_enrichment", "H transition enrichment", CONTEXT, "expression-score difference", False, "", "hv", 2, True),
    MetricDefinition("H_spatial_variance", "H spatial variance", CONTEXT, "variance", False, "", "hv", 3, True),
    MetricDefinition("H_coefficient_of_variation", "H coefficient of variation", CONTEXT, "ratio", False, "", "hv", 4, True),
    MetricDefinition("H_local_hotspot_fraction", "H local-hotspot fraction", CONTEXT, "fraction", False, "pooled_reference_target_q90", "hv", 3, True),
    MetricDefinition("V_raw_mean", "V raw-scale mean", CONTEXT, "expression score", False, "", "hv", 3, True),
    MetricDefinition("V_raw_median", "V raw-scale median", CONTEXT, "expression score", False, "", "hv", 2, True),
    MetricDefinition("V_q75", "V 75th percentile", CONTEXT, "expression score", False, "", "hv", 3, True),
    MetricDefinition("V_q90", "V 90th percentile", CONTEXT, "expression score", False, "", "hv", 2, True),
    MetricDefinition("V_high_fraction", "V pooled-threshold high fraction", CONTEXT, "fraction", False, "pooled_reference_target_q90", "hv", 2, True),
    MetricDefinition("V_variance", "V variance", CONTEXT, "variance", False, "", "hv", 3, True),
    MetricDefinition("V_MAD", "V median absolute deviation", CONTEXT, "expression score", False, "", "hv", 3, True),
    MetricDefinition("V_transition_enrichment", "V transition enrichment", CONTEXT, "expression-score difference", False, "", "hv", 2, True),
    MetricDefinition("V_spatial_variance", "V spatial variance", CONTEXT, "variance", False, "", "hv", 3, True),
    MetricDefinition("V_coefficient_of_variation", "V coefficient of variation", CONTEXT, "ratio", False, "", "hv", 4, True),
    MetricDefinition("V_local_hotspot_fraction", "V local-hotspot fraction", CONTEXT, "fraction", False, "pooled_reference_target_q90", "hv", 3, True),
)


REGISTRY_BY_NAME = {definition.internal_name: definition for definition in METRIC_REGISTRY}
DELTA_METRIC_SPECS = tuple(
    (definition.resolved_delta_name, definition.internal_name) for definition in METRIC_REGISTRY
)
GROUP_METRICS = tuple(definition.internal_name for definition in METRIC_REGISTRY)


def metric_definition(name: str) -> MetricDefinition:
    return REGISTRY_BY_NAME.get(
        name,
        MetricDefinition(name, name.replace("_", " ").title(), "Other", "value"),
    )


def metrics_for_plot_group(plot_group: str) -> tuple[MetricDefinition, ...]:
    return tuple(definition for definition in METRIC_REGISTRY if definition.plot_group == plot_group)


def heatmap_metrics() -> tuple[MetricDefinition, ...]:
    return tuple(
        definition
        for definition in METRIC_REGISTRY
        if not definition.deprecated
        and definition.interpretation_priority <= 3
        and definition.plot_group not in {"sample_scale", "topology_raw"}
    )


def registry_dataframe_records() -> list[dict]:
    return [definition.to_dict() for definition in METRIC_REGISTRY]
