"""Example of using an imported geometry in a PorePy model.

Loads fracture network csv file `shared/disks.csv` and runs a single-phase
flow model on it.

"""

from pathlib import Path

import porepy as pp
from porepy.models.fluid_mass_balance import SinglePhaseFlow


class ImportedGeometry:
    # params: dict
    units: pp.Units

    def set_domain(self) -> None:
        """Set the domain based on the CSV."""
        self.create_fracture_network()
        self._domain = self.fracture_network.domain

    def set_fractures(self) -> None:
        self.create_fracture_network()
        self._fractures = self.fracture_network.fractures

    def create_fracture_network(self) -> None:
        """Set the fracture network from the CSV geometry file."""
        csv_geometry_file = Path("shared/disks.csv")
        expected_domain_size = 1000
        # TODO: read domain size from csv file (first line)

        self.fracture_network = pp.fracture_importer.network_from_csv(
            csv_geometry_file,
            has_domain=True,
            tol=expected_domain_size * 1e-6,
        )

    def grid_type(self) -> str:
        return "simplex"

    def meshing_arguments(self) -> dict:
        mesh_args = {}
        mesh_args["cell_size"] = self.units.convert_units(100.0, "m")
        mesh_args["cell_size_fracture"] = self.units.convert_units(100.0, "m")
        return mesh_args


class SinglePhaseFlowGeometry(ImportedGeometry, SinglePhaseFlow):
    """Combining the imported geometry and the default model."""


model_params = {}
model = SinglePhaseFlowGeometry(model_params)
pp.ModelRunner(model).run()
