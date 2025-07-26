# Genomics Project - Path Configuration

This project uses relative paths to ensure portability across different systems and users.

## Directory Structure

When you clone this repository, the expected structure is:

```
your-project-folder/
├── web_app/
│   ├── shiny_ui.py
│   ├── utilities.py
│   ├── custom_css.py
│   ├── old_intro_consent_card_ui.py
│   └── www/
│       ├── blood5.jpg
│       ├── sports2.jpg
│       ├── gene2.webp
│       ├── pexels-pixabay-531880.jpg
│       └── other_images...
├── background_data/
│   └── sports_snps.txt
├── genomes/
│   ├── genome_Cajun_v5_Full_20231121192441.txt
│   ├── genome_Frederik_FangelTolberg_v5_Full_20241117223640.txt
│   ├── genome_Marika_Forsythe_v4_Full_20240828221950.txt
│   └── hu278AF5_20210124151934.txt
└── gwas_data_filtered_phased_low_filt.csv
```

## Path Resolution

The application now uses relative paths based on the script location:

- **Images**: Located relative to `web_app/www/`
- **GWAS Data**: Located relative to project root
- **Background Data**: Located relative to project root in `background_data/`
- **Genome Files**: Located relative to project root in `genomes/`

## How It Works

1. **Script Directory Detection**: The code automatically detects where the Python files are located using `Path(__file__).parent`

2. **Relative Path Construction**: All file paths are constructed relative to the detected script location

3. **Cross-Platform Compatibility**: Uses `pathlib.Path` which works on Windows, Mac, and Linux

## For Developers

If you need to add new file references:

```python
from pathlib import Path

# Get current script directory
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent

# Reference files relative to project structure
image_path = SCRIPT_DIR / "www" / "your_image.jpg"
data_path = PROJECT_ROOT / "your_data.csv"
genome_path = PROJECT_ROOT / "genomes" / "your_genome.txt"
```

## Helper Functions

The `utilities.py` file now includes helper functions:

- `get_genome_file_path(filename)`: Get full path to a genome file
- `list_available_genome_files()`: List all available genome files

## Environment Variables (Optional)

You can also set environment variables to override default paths if needed:

```bash
export GENOMICS_PROJECT_ROOT=/path/to/your/project
export GENOMICS_GENOMES_DIR=/path/to/genomes
export GENOMICS_GWAS_DATA=/path/to/gwas_data.csv
```

The application will use these if set, otherwise it defaults to relative paths.
