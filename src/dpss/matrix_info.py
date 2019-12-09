from typing import (
    Optional,
    FrozenSet,
)

from attr import dataclass
from more_itertools import one


@dataclass
class MatrixInfo:
    source: str
    project_uuid: str
    zip_path: Optional[str]
    extract_path: str
    lib_con_approaches: FrozenSet[str] = frozenset()

    @property
    def figures_folder(self) -> str:
        lca_count = len(self.lib_con_approaches)
        if lca_count == 1:
            from dpss.matrix_summary_stats import MatrixSummaryStats
            lca = one(self.lib_con_approaches)
            suffix = MatrixSummaryStats.translate_lca(lca)
            return f'{self.project_uuid}/{suffix}/'
        else:
            raise RuntimeError(f'Should not upload figures for matrix {self.project_uuid}'
                               'because it has not been separated by library construction approach.')
