#!/usr/bin/env python3
"""
Minimal Jupyter to Quarto Converter
"""

import json
import os
from pathlib import Path

def convert_notebooks_simple():
    """Simple conversion without external dependencies"""
    
    print("Looking for Jupyter notebooks...")
    
    # Find all notebooks
    notebooks = list(Path('presentation_files').glob('*.ipynb'))
    
    if not notebooks:
        print("No .ipynb files found in current directory")
        return
    
    presentation = """---
title: "Modflow 6 Lab: Customizing for Flow, Transport and Geothermal Modeling"
format: revealjs
---

# Workshop 

"""
    
    for notebook_path in notebooks:
        print(f"Converting: {notebook_path}")
        
        try:
            with open(notebook_path, 'r', encoding='utf-8') as f:
                nb = json.load(f)
            
            presentation += f"## {notebook_path.stem.replace('_', ' ').title()}\n\n"
            
            for cell in nb['cells']:
                if cell['cell_type'] == 'markdown':
                    presentation += ''.join(cell['source']) + '\n\n'
                elif cell['cell_type'] == 'code':
                    code = ''.join(cell['source']).strip()
                    if code:
                        presentation += f"```python\n{code}\n```\n\n"
                        
        except Exception as e:
            print(f"Error processing {notebook_path}: {e}")
    
    # Add conclusions
    presentation += """
## Conclusions

.....
---
## Thank You
"""
    
    with open('presentation.qmd', 'w', encoding='utf-8') as f:
        f.write(presentation)
    
    print("✅ Created presentation.qmd")
    print("Next: Install Quarto and run: quarto preview presentation.qmd")

if __name__ == "__main__":
    convert_notebooks_simple()

