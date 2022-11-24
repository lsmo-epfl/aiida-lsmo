from pathlib import Path
from typing import Optional

from aiida_testing.mock_code import InputHasher


class CustomInputHasher(InputHasher):

    def __call__(self, cwd: Path) -> str:
        hash_string = super().__call__(cwd)
        return hash_string
