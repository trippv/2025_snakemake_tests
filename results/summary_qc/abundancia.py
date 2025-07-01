from multiqc.utils import config
from multiqc.plots import heatmap
from multiqc.modules.base_module import BaseMultiqcModule

import random

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Inicializa el módulo con su nombre y configuración
        super().__init__(
            name="Abundancia",
            anchor="abundancia",
            target="Expresión génica simulada",
            href="",
            info="Este módulo muestra un heatmap de datos de expresión génica simulados."
        )

        # Simula los datos
        genes = [f"Gene_{i}" for i in range(1, 11)]
        samples = [f"Sample_{j}" for j in range(1, 7)]
        heatmap_data = {
            gene: {sample: round(random.uniform(0, 100), 2) for sample in samples}
            for gene in genes
        }

        # Configura el heatmap
        pconfig = {
            "id": "abundancia_heatmap",
            "title": "Abundancia: Heatmap de expresión génica",
            "xlab": "Muestras",
            "ylab": "Genes",
            "tt_label": "{x}, {y}: {value:.2f}",
            "height": 500
        }

        # Añade la sección con el heatmap
        self.add_section(
            name="Expresión génica",
            anchor="expresion_genica",
            description="Este heatmap muestra la abundancia relativa de genes simulados en diferentes muestras.",
            helptext="Valores generados aleatoriamente para representar niveles de expresión.",
            plot=heatmap.plot(data=heatmap_data, xcats=samples, ycats=genes, pconfig=pconfig)
        )
