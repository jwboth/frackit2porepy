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
        self._fractures = (
            self.fracture_network.fractures # + self.quad_fracture_network.fractures
        )

    def create_fracture_network(self) -> None:
        """Set the fracture network from the CSV geometry file."""
        self.fracture_network = pp.fracture_importer.network_from_csv(
            Path("shared/disks.csv"),
            has_domain=True,
            tol=1e-3,  # NOTE: Small wrt characteristic fracture/domain size.
        )
        # self.quad_fracture_network = pp.fracture_importer.network_from_csv(
        #     Path("shared/quads.csv"),
        #     has_domain=False,
        #     tol=1000,  # NOTE: Small wrt characteristic fracture/domain size.
        # )

    def grid_type(self) -> str:
        """Use simplex meshes to handle the complex geometry."""
        return "simplex"

    def meshing_arguments(self) -> dict:
        """Assign coarse mesh size for the imported geometry."""
        mesh_args = {}
        mesh_args["cell_size"] = self.units.convert_units(100.0, "m")
        mesh_args["cell_size_fracture"] = self.units.convert_units(100.0, "m")
        return mesh_args


class SinglePhaseFlowGeometry(ImportedGeometry, SinglePhaseFlow):
    """Combining the imported geometry and the default model."""


model_params = {}
model = SinglePhaseFlowGeometry(model_params)
pp.ModelRunner(model).run()
