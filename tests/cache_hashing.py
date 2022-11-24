from pathlib import Path
from typing import Optional

from aiida_testing.mock_code import InputHasher


class CustomInputHasher(InputHasher):

    def __call__(self, cwd: Path) -> str:
        hash_string = super().__call__(cwd)
        if hash_string == '805a68cedb2a2d30bab0ee403bc9abed':
            # temporary fix for: https://github.com/lsmo-epfl/aiida-lsmo/issues/102
            return 'e323f1f24dd00a44b9c973208dcd63d8'
        return hash_string
