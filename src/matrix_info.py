from typing import Optional

from attr import dataclass


@dataclass
class MatrixInfo:
    source: str
    project_uuid: str
    zip_path: Optional[str]
    extract_path: str
    lib_prep_method: Optional[str] = None
