#%%
import importlib
from pathlib import Path
import os
    
for path in Path(__file__).parent.glob('*.py'):
    if path.name == "__init__.py":
        continue
    importlib.import_module(f".{path.stem}", package = __package__)