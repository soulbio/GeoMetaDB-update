# GEOmetaDB Database Update Tool

This tool updates the GEOmetaDB SQLite database with missing GSE IDs and their comprehensive metadata.

## Features

The tool updates three main tables in the GEOmetaDB SQLite database:

### 1. `gse` table
Updates the following metadata for each GSE series:
- `title`: Series title
- `summary`: Series summary/description  
- `overall_design`: Overall experimental design
- `pubmed_id`: Associated PubMed ID
- `type`: Series type

### 2. `gsm` table  
Updates the following metadata for each GSM sample:
- `title`: Sample title
- `source_name_ch1`: Source name
- `gpl`: Platform ID
- `organism_ch1`: Organism information
- `characteristics_ch1`: Sample characteristics

### 3. `gse_gsm` table
Creates relationships between GSE series and their GSM samples:
- `gse`: GSE series ID
- `gsm`: GSM sample ID

## Usage

### Command Line Usage

```bash
python GEOmetaDB_update.py --sqlite_path /path/to/GeoMetaDB.sqlite 
```

#### Optional Parameters:
- `--batch_size 10`: Number of GSE IDs to process in each batch (default: 10)
- `--delay 1.0`: Delay in seconds between API calls to avoid rate limiting (default: 1.0)

### Example Usage:

```bash
# Basic usage
python GEOmetaDB_update.py --sqlite_path ./GeoMetaDB.sqlite 

# With custom batch size and delay
python GEOmetaDB_update.py --sqlite_path ./GeoMetaDB.sqlite --batch_size 5 --delay 2.0
```